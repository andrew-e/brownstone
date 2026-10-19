fixture <- function(f) readLines(test_path("fixtures", f))
url <- "https://ldlink.nih.gov/LDlinkRest/ldproxy?var=rs1&pop=EUR&token=SECRET123"

# What base R's url() reports for an HTTP error: a warning with the status and
# an error, both quoting the URL, token included.
http_failure <- function(status) {
  function(u) {
    warning("cannot open URL '", u, "': HTTP status was '", status, "'")
    stop("cannot open the connection to '", u, "'")
  }
}

test_that("an LDproxy response is parsed into standard columns", {
  ld <- .parse_ldproxy(fixture("ldproxy_rs11591147_20kb.txt"))
  expect_named(ld, c("chr", "pos", "snp", "alleles", "maf", "distance", "dprime", "r2"))
  q <- ld[ld$distance == 0, ]
  expect_equal(q$snp, "rs11591147")
  expect_equal(q$chr, "1")
  expect_equal(q$pos, 55039974)
  expect_equal(q$r2, 1)
  expect_true(all(ld$r2 >= 0.01))           # LDproxy leaves out r2 <= 0.01 (rounded)
  expect_true(all(abs(ld$distance) <= 2e4)) # the window is +/- around the variant
})

test_that("LDlink errors become R errors", {
  expect_null(.ldlink_error(fixture("ldproxy_rs11591147_20kb.txt")))
  expect_error(.parse_ldproxy(fixture("ldproxy_error.txt")), "Invalid or expired API token")
  expect_error(.parse_ldproxy("RS_Number\tR2"), "unexpected table")
})

test_that("variants are given to LDlink as rsIDs or chr:pos", {
  expect_equal(.ldlink_var("RS7412"), "rs7412")
  expect_equal(.ldlink_var("19:44908822"), "chr19:44908822")
  expect_equal(.ldlink_var("19:44908822:C:T"), "chr19:44908822")
  expect_equal(.ldlink_var("chr23:100"), "chrX:100")
  expect_error(.ldlink_var("Y:100"), "no LD")
  expect_error(.ldlink_var("APOE"), "rsID")
})

test_that("a busy or failing service is retried with backoff, and the token never shows", {
  calls <- 0
  fail <- http_failure("429 Too Many Requests")
  local_mocked_bindings(.read_url = function(u) {
    calls <<- calls + 1
    if (calls < 3) fail(u) else fixture("ldproxy_rs11591147_20kb.txt")
  })
  msgs <- testthat::capture_messages(out <- .ldlink_get(url, wait = 0))
  expect_equal(calls, 3)
  expect_length(msgs, 2)
  expect_match(msgs[1], "429 Too Many Requests")
  expect_false(any(grepl("SECRET123", msgs)))
  expect_true(startsWith(out[1], "RS_Number"))

  # A "busy" message from LDlink itself is retried too.
  calls <- 0
  local_mocked_bindings(.read_url = function(u) {
    calls <<- calls + 1
    if (calls == 1) '{"error": "Concurrent API requests restricted for this token"}' else
      fixture("ldproxy_rs11591147_20kb.txt")
  })
  expect_message(.ldlink_get(url, wait = 0), "Concurrent")
  expect_equal(calls, 2)
})

test_that("other failures stop at once, and retries give up", {
  calls <- 0
  fail <- http_failure("404 Not Found")
  local_mocked_bindings(.read_url = function(u) {
    calls <<- calls + 1
    fail(u)
  })
  err <- expect_error(.ldlink_get(url, wait = 0), "404 Not Found")
  expect_equal(calls, 1)
  expect_false(grepl("SECRET123", conditionMessage(err)))

  calls <- 0
  fail <- http_failure("503 Service Unavailable")
  local_mocked_bindings(.read_url = function(u) {
    calls <<- calls + 1
    fail(u)
  })
  err <- expect_error(suppressMessages(.ldlink_get(url, tries = 3, wait = 0)),
                      "did not answer after 3 tries")
  expect_equal(calls, 3)
  expect_false(grepl("SECRET123", conditionMessage(err)))

  local_mocked_bindings(.read_url = function(u) fixture("ldproxy_error.txt"))
  expect_error(.ldlink_get(url, wait = 0), "Invalid or expired API token")
})

test_that("ld_proxy() asks LDlink once per variant, population and build", {
  rm(list = ls(.ld_cache), envir = .ld_cache)
  urls <- character(0)
  local_mocked_bindings(.ldlink_get = function(u, ...) {
    urls <<- c(urls, u)
    fixture("ldproxy_rs11591147_20kb.txt")
  })
  a <- ld_proxy("rs11591147", window = 2e4, token = "SECRET123")
  b <- ld_proxy("rs11591147", window = 1e4, token = "SECRET123")   # from the cache
  expect_length(urls, 1)
  expect_true(all(abs(b$distance) <= 1e4))
  expect_lt(nrow(b), nrow(a))
  expect_match(urls[1], "var=rs11591147&pop=EUR&r2_d=r2&window=20000", fixed = TRUE)
  expect_match(urls[1], "genome_build=grch38_high_coverage", fixed = TRUE)
  ld_proxy("rs11591147", window = 2e4, pop = "FIN", token = "SECRET123")
  ld_proxy("rs11591147", window = 2e4, build = "GRCh37", token = "SECRET123")
  expect_length(urls, 3)
  expect_match(urls[3], "genome_build=grch37", fixed = TRUE)

  expect_error(ld_proxy("rs1", token = ""), "access token")
  expect_error(ld_proxy("rs1", window = 2e6, token = "x"), "window")
  rm(list = ls(.ld_cache), envir = .ld_cache)
})

test_that("GRCh38 falls back to the original release when the high-coverage one fails", {
  rm(list = ls(.ld_cache), envir = .ld_cache)
  urls <- character(0)
  local_mocked_bindings(.ldlink_get = function(u, ...) {
    urls <<- c(urls, u)
    if (grepl("high_coverage", u)) '{"error": "Variant is monoallelic in the chosen population(s)."}'
    else fixture("ldproxy_rs11591147_20kb.txt")
  })
  expect_message(ld <- ld_proxy("rs11591147", window = 2e4, token = "x"), "original 1000 Genomes")
  expect_length(urls, 2)
  expect_match(urls[2], "genome_build=grch38&", fixed = TRUE)
  expect_gt(nrow(ld), 0)
  # Unreachable LDlink is not retried against another release.
  rm(list = ls(.ld_cache), envir = .ld_cache)
  urls <- character(0)
  local_mocked_bindings(.ldlink_get = function(u, ...) {
    urls <<- c(urls, u)
    stop("Could not reach LDlink (HTTP '404 Not Found').", call. = FALSE)
  })
  expect_error(ld_proxy("rs11591147", window = 2e4, token = "x"), "Could not reach")
  expect_length(urls, 1)
  rm(list = ls(.ld_cache), envir = .ld_cache)
})

test_that("LD rows are matched by ID, then by an unambiguous position", {
  d <- standardise_sumstats(data.frame(CHR = c(1, 1, 1, 23), BP = c(100, 200, 200, 300),
                                       SNP = c("rs1", "a", "b", "rsX"), P = 0.5))
  tab <- data.frame(snp = c("rs1", ".", NA, "rs9"), chr = c("1", "1", "1", "X"),
                    pos = c(999, 200, 100, 300))
  expect_equal(.match_variants(tab, d, unique_pos = TRUE), c(1L, NA, 1L, 4L))
  # Labels keep their first match at a shared position.
  expect_equal(.match_variants(tab, d), c(1L, 2L, 1L, 4L))
})

test_that("LDlink answers (live)", {
  skip_on_cran()
  skip_if_offline("ldlink.nih.gov")
  skip_if_not(nzchar(Sys.getenv("LDLINK_TOKEN")), "no LDLINK_TOKEN")
  ld <- ld_proxy("rs7412", window = 1e4)
  expect_equal(ld$r2[ld$distance == 0], 1)
  expect_equal(ld$pos[ld$distance == 0], 44908822)
})
