# Five independent loci, one variant each, plus noise below every line: three
# genuine effects whose case/control frequencies agree with the reported
# beta, one whose estimate its own allele counts don't support (the MAN2B1
# case), and one whose crude odds ratio disagrees in direction but by less
# than either flag margin.
artifact_gwas <- function() {
  f <- c(0.30, 0.20, 0.10, 0.01, 0.05)
  beta <- c(0.8, 1.5, -0.6, 12.6, -0.165)
  se <- c(0.08, 0.10, 0.09, 1.61, 0.03)
  or <- c(exp(beta[1:3]), 1.75, 1.01)
  p2 <- f
  p1 <- or * p2 / (1 - p2 + or * p2)
  n <- 20
  data.frame(
    CHR = c(as.character(1:5), rep("6", n)),
    BP = c(seq_along(f) * 1e6 + 5e5, seq_len(n) * 1e6),
    P = c(signif(2 * stats::pnorm(-abs(beta / se)), 4), rep(0.5, n)),
    BETA = c(beta, rep(0, n)), SE = c(se, rep(0.1, n)),
    EAF = c(f, rep(0.25, n)),
    af_case = c(p1, rep(0.25, n)), af_ctrl = c(p2, rep(0.25, n))
  )
}

n_ring <- function(p) {
  sum(vapply(p$layers, function(l) {
    inherits(l$geom, "GeomPoint") && identical(l$aes_params$shape, 1) &&
      identical(l$aes_params$colour, .ink$primary)
  }, logical(1)))
}

test_that("each lead is compared with its own case/control frequencies", {
  d <- standardise_sumstats(artifact_gwas(), ccafs = TRUE)
  out <- .artifact_table(d, 5e-8, 5e5)
  expect_length(out$flagged, 5)
  i <- match("4:4500000", out$snp)
  expect_equal(out$crude_or[i], 1.75)
  expect_equal(out$discrepancy[i], 12.6 - log(1.75))
  expect_equal(out$z_discrepancy[i], (12.6 - log(1.75)) / 1.61)
  expect_equal(out$ratio[i], exp(12.6 - log(1.75)))
  # the genuine effects: their frequencies imply the reported estimate
  genuine <- match(c("1:1500000", "2:2500000", "3:3500000"), out$snp)
  expect_true(all(abs(out$discrepancy[genuine]) < 1e-12))
  # a direction flip alone is not flagged: with few carriers the crude odds
  # ratio is noise
  j <- match("5:5500000", out$snp)
  expect_true(out$direction_disagrees[j])
  expect_equal(out$z_discrepancy[j], (0.165 + log(1.01)) / 0.03)
  expect_false(out$flagged[j])
  # only the estimate its own counts don't support
  expect_equal(which(out$flagged), i)
})

test_that("find_artifacts returns the numbers, the genes and the message", {
  expect_message(out <- find_artifacts(artifact_gwas()), "see find_artifacts")
  expect_s3_class(out, "data.frame")
  expect_false("row" %in% names(out))
  expect_true("gene" %in% names(out))
  expect_equal(sum(out$flagged), 1)
  # unannotated
  out <- suppressMessages(find_artifacts(artifact_gwas(), annotate = NULL))
  expect_false("gene" %in% names(out))
})

test_that("frequencies of 0, 1 or missing cannot be judged", {
  x <- artifact_gwas()
  x$af_case[1] <- 0   # may only be rounding
  x$af_ctrl[2] <- 1
  x$af_case[3] <- NA
  out <- suppressMessages(find_artifacts(x))
  expect_equal(nrow(out), 2)
})

test_that("frequencies are folded to the same orientation as beta", {
  f <- 0.3; r <- 2
  p2m <- f
  p1m <- r * p2m / (1 - p2m + r * p2m)
  # reported against the major allele: beta is flipped to the minor allele
  # and the frequencies with it, so the crude odds ratio still agrees
  x <- data.frame(CHR = "7", BP = 1e6, P = 1e-30,
                  BETA = -log(r), SE = 0.01, EAF = 1 - f,
                  af_case = 1 - p1m, af_ctrl = 1 - p2m)
  d <- standardise_sumstats(x, ccafs = TRUE)
  expect_equal(d$beta, log(r))
  expect_equal(d$af_case, p1m)
  expect_equal(d$af_ctrl, p2m)
  out <- .artifact_table(d, 5e-8, 5e5)
  expect_length(out$flagged, 1)
  expect_true(abs(out$discrepancy) < 1e-12)
  expect_false(out$flagged)
})

test_that("the margins are the arguments", {
  d <- standardise_sumstats(artifact_gwas(), ccafs = TRUE)
  # the direction flip is 5.8 SEs: a z of 5 flags it
  expect_equal(sum(.artifact_table(d, 5e-8, 5e5, z = 5)$flagged), 2)
  expect_equal(sum(.artifact_table(d, 5e-8, 5e5, z = Inf)$flagged), 1)
  expect_equal(sum(.artifact_table(d, 5e-8, 5e5, ratio = Inf)$flagged), 1)
})

test_that("the plots mark flagged leads and stay quiet without frequencies", {
  x <- artifact_gwas()
  pb <- suppressMessages(brownstone(x, labels = NULL))
  expect_equal(n_ring(pb[[1]]), 0)   # the main panel: no marker
  expect_equal(n_ring(pb[[2]]), 1)   # the effect track
  expect_message(p <- effect_vs_maf(x, labels = NULL), "see find_artifacts")
  expect_equal(n_ring(p), 1)
  expect_message(p <- volcano_plot(x, labels = NULL), "see find_artifacts")
  expect_equal(n_ring(p), 1)
  # no case/control frequency columns: no marker, no message, no error
  gw <- simulate_gwas(n_variants = 5e4, seed = 1)
  pb <- expect_silent(brownstone(gw, labels = NULL))
  expect_equal(n_ring(pb[[2]]), 0)
  expect_equal(n_ring(effect_vs_maf(gw, labels = NULL)), 0)
  expect_silent(volcano_plot(gw, labels = NULL))
})

test_that("standardise carries case/control frequencies only when asked", {
  x <- artifact_gwas()
  d <- standardise_sumstats(x)
  expect_false("af_case" %in% names(d))
  d <- standardise_sumstats(x, ccafs = TRUE)
  expect_true(all(c("af_case", "af_ctrl") %in% names(d)))
  # SAIGE's names (AF_case/AF_ctrl) are found case-insensitively
  names(x)[match(c("af_case", "af_ctrl"), names(x))] <- c("AF_case", "AF_ctrl")
  d <- standardise_sumstats(x, ccafs = TRUE)
  expect_true(all(c("af_case", "af_ctrl") %in% names(d)))
  # FinnGen release 12's names
  names(x)[match(c("AF_case", "AF_ctrl"), names(x))] <-
    c("af_alt_cases", "af_alt_controls")
  d <- standardise_sumstats(x, ccafs = TRUE)
  expect_true(all(c("af_case", "af_ctrl") %in% names(d)))
  # one frequency alone is not carried
  x$af_alt_controls <- NULL
  d <- standardise_sumstats(x, ccafs = TRUE)
  expect_false("af_case" %in% names(d))
})
