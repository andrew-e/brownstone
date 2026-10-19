#' Linkage disequilibrium with one variant, from LDlink
#'
#' @description
#' Asks the LDproxy service of LDlink (<https://ldlink.nih.gov>) for the
#' linkage disequilibrium (LD) between one variant and every variant within
#' `window` base pairs of it, in a 1000 Genomes reference population.
#' [region_plot()] calls it for you; call it yourself to fetch LD once and
#' pass it to several plots.
#'
#' * **A personal access token is needed.** Register (free) on the API Access
#'   tab of LDlink and add `LDLINK_TOKEN=<token>` to your `.Renviron`. The token
#'   never appears in messages or errors.
#' * **Only variants with r² > 0.01 are returned** by LDproxy. A variant that is
#'   missing is either in negligible LD with `snp` or not in the 1000 Genomes
#'   reference.
#' * **The reference** for GRCh38 is the high-coverage (30x) 1000 Genomes
#'   release, which has more rare variants than the original; if LDlink can't
#'   use the variant there (it calls some common variants monoallelic), the
#'   original release. For GRCh37, the original release.
#' * **Results are cached** for the R session, so plotting the same region
#'   again doesn't call LDlink again.
#' * **A busy or failing service** (HTTP 429 or 5xx, timeouts) is retried with
#'   exponential backoff: after 2, 4 and 8 seconds.
#'
#' @param snp The variant: an rsID (`"rs7412"`), or a position (`"19:44908822"`
#'   or `"chr19:44908822"`) in `build` coordinates.
#' @param pop 1000 Genomes population code(s), e.g. `"EUR"` (the default),
#'   `"FIN"` or `c("CEU", "GBR")`.
#' @param build Genome build of the positions, `"GRCh38"` or `"GRCh37"`.
#' @param window Base pairs either side of `snp`, at most 1,000,000.
#' @param token Your LDlink access token.
#'
#' @return A data frame with one row per variant in LD with `snp`: `chr`,
#'   `pos`, `snp` (rsID, `NA` if none), `alleles`, `maf` (in the reference
#'   population), `distance` (base pairs from `snp`), `dprime` and `r2`.
#' @export
#' @examples
#' if (nzchar(Sys.getenv("LDLINK_TOKEN"))) {
#'   ld <- ld_proxy("rs7412", window = 1e5)
#'   head(ld)
#' }
ld_proxy <- function(snp, pop = "EUR", build = c("GRCh38", "GRCh37"), window = 5e5,
                     token = Sys.getenv("LDLINK_TOKEN")) {
  build <- match.arg(build)
  if (!nzchar(token)) {
    stop("LDlink needs a personal access token: register on the API Access tab of ",
         "https://ldlink.nih.gov and add LDLINK_TOKEN=<token> to your .Renviron.",
         call. = FALSE)
  }
  if (!is.numeric(window) || length(window) != 1L || window < 0 || window > 1e6) {
    stop("`window` must be between 0 and 1,000,000 base pairs.", call. = FALSE)
  }
  var <- .ldlink_var(snp)
  window <- as.integer(ceiling(window))
  pop <- toupper(pop)
  key <- paste(var, paste(pop, collapse = "+"), build, sep = "|")
  hit <- .ld_cache[[key]]
  if (!is.null(hit) && hit$window >= window) {
    return(hit$ld[abs(hit$ld$distance) <= window, , drop = FALSE])
  }
  # For GRCh38, the high-coverage release first; LDlink's own errors there
  # (e.g. a variant it calls monoallelic) fall back to the original release.
  refs <- if (build == "GRCh38") c("grch38_high_coverage", "grch38") else "grch37"
  for (k in seq_along(refs)) {
    url <- sprintf("%s/ldproxy?var=%s&pop=%s&r2_d=r2&window=%d&genome_build=%s&token=%s",
                   .ldlink_root, var, paste(pop, collapse = "%2B"), window, refs[k], token)
    ld <- tryCatch(.parse_ldproxy(.ldlink_get(url)), ldlink_error = function(e) {
      if (k == length(refs)) stop(e)
      message(conditionMessage(e), " Trying the original 1000 Genomes GRCh38 release.")
      NULL
    })
    if (!is.null(ld)) break
  }
  assign(key, list(window = window, ld = ld), envir = .ld_cache)
  ld
}

# An error reported by LDlink itself (a bad token, an unknown variant), as
# opposed to LDlink being unreachable.
.ldlink_stop <- function(msg) {
  stop(structure(class = c("ldlink_error", "error", "condition"),
                 list(message = paste0("LDlink: ", .hide_token(msg)), call = NULL)))
}

.ldlink_root <- "https://ldlink.nih.gov/LDlinkRest"

# LDproxy results for this R session, by variant, population and build.
.ld_cache <- new.env(parent = emptyenv())

# The variant as LDlink wants it: an rsID, or "chrN:pos".
.ldlink_var <- function(snp) {
  snp <- trimws(as.character(snp))
  if (length(snp) != 1L || is.na(snp)) stop("`snp` must be one variant.", call. = FALSE)
  if (grepl("^rs[0-9]+$", snp, ignore.case = TRUE)) return(tolower(snp))
  m <- regmatches(snp, regexec("^(chr)?([0-9]{1,2}|X|Y|XY|MT|M)[:_]([0-9]+)", snp,
                               ignore.case = TRUE))[[1]]
  if (!length(m)) {
    stop("`snp` must be an rsID (\"rs7412\") or a position (\"19:44908822\"), not \"",
         snp, "\".", call. = FALSE)
  }
  chr <- .norm_chr(m[3])
  if (!chr %in% c(as.character(1:22), "X")) {
    stop("LDlink has no LD for chromosome ", chr, ".", call. = FALSE)
  }
  paste0("chr", chr, ":", m[4])
}

.hide_token <- function(x) gsub("token=[^&'\" ]*", "token=<hidden>", x)

# Read one LDlink URL. Retries with exponential backoff (2, 4, 8 s) when the
# service is busy or failing: HTTP 429 or 5xx, a timeout, or an error message
# saying so. Other failures (a bad token, an unknown variant) stop at once.
# Messages never include the token.
.ldlink_get <- function(url, tries = 4, wait = 2) {
  old <- options(timeout = max(300, getOption("timeout")))
  on.exit(options(old), add = TRUE)
  problem <- "no response"
  for (i in seq_len(tries)) {
    warnings <- character(0)
    out <- withCallingHandlers(
      tryCatch(.read_url(url), error = function(e) structure(conditionMessage(e), class = "ldlink_fail")),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    if (!inherits(out, "ldlink_fail")) {
      err <- .ldlink_error(out)
      if (is.null(err)) return(out)
      if (!grepl("busy|concurrent|too many|try again|timed? ?out", err, ignore.case = TRUE)) {
        .ldlink_stop(err)
      }
      problem <- err
    } else {
      txt <- .hide_token(c(warnings, unclass(out)))
      status <- unlist(regmatches(txt, regexpr("HTTP status was '[^']*'", txt)))
      code <- suppressWarnings(as.integer(substr(sub("HTTP status was '", "", status), 1, 3)))
      transient <- any(code == 429 | code >= 500, na.rm = TRUE) ||
        any(grepl("Timeout was reached|timed out", txt, ignore.case = TRUE))
      problem <- if (length(status)) sub("HTTP status was ", "HTTP ", status[1]) else txt[length(txt)]
      if (!transient) stop("Could not reach LDlink (", problem, ").", call. = FALSE)
    }
    if (i < tries) {
      pause <- wait * 2^(i - 1)
      message("LDlink: ", .hide_token(problem), "; trying again in ", pause, " s.")
      Sys.sleep(pause)
    }
  }
  stop("LDlink did not answer after ", tries, " tries (", .hide_token(problem), ").",
       call. = FALSE)
}

.read_url <- function(url) {
  con <- url(url)
  on.exit(close(con))
  readLines(con, warn = FALSE, encoding = "UTF-8")
}

# LDlink reports errors as JSON ({"error": "..."}) or as text in place of the
# table. NULL when `lines` is a table.
.ldlink_error <- function(lines) {
  if (!length(lines)) return("an empty response")
  if (startsWith(lines[1], "RS_Number")) return(NULL)
  txt <- paste(lines, collapse = " ")
  m <- regmatches(txt, regexec("\"error\"\\s*:\\s*\"([^\"]*)\"", txt))[[1]]
  if (length(m)) return(m[2])
  if (grepl("error", txt, ignore.case = TRUE)) return(trimws(substr(txt, 1, 300)))
  NULL
}

# An LDproxy table as a data frame with standard names.
.parse_ldproxy <- function(lines) {
  err <- .ldlink_error(lines)
  if (!is.null(err)) .ldlink_stop(err)
  d <- utils::read.delim(text = lines, stringsAsFactors = FALSE, quote = "",
                         na.strings = c("NA", "."), check.names = FALSE)
  if (!all(c("Coord", "R2") %in% names(d))) {
    stop("LDlink returned an unexpected table (columns: ", paste(names(d), collapse = ", "),
         ").", call. = FALSE)
  }
  coord <- sub("^chr", "", d$Coord, ignore.case = TRUE)
  col <- function(nm, na) if (nm %in% names(d)) d[[nm]] else rep(na, nrow(d))
  data.frame(chr = .norm_chr(sub(":.*$", "", coord)),
             pos = as.numeric(sub("^.*:", "", coord)),
             snp = as.character(col("RS_Number", NA_character_)),
             alleles = as.character(col("Alleles", NA_character_)),
             maf = as.numeric(col("MAF", NA_real_)),
             distance = as.numeric(col("Distance", NA_real_)),
             dprime = as.numeric(col("Dprime", NA_real_)),
             r2 = as.numeric(d$R2), stringsAsFactors = FALSE)
}
