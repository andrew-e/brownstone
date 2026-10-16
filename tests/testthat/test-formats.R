# Output formats of the common GWAS tools. Headers are copied from each tool's
# documented output; the values are simulated.

make_rows <- function(n = 6) {
  set.seed(1)
  data.frame(chr = rep(c(1, 23), length.out = n), pos = seq_len(n) * 1000,
             id = paste0("rs", seq_len(n)), af = c(0.2, 0.7, 0.001, 0.4, 0.05, 0.9)[seq_len(n)],
             beta = c(0.3, -0.2, 1.5, 0.1, -0.4, 0.2)[seq_len(n)], se = 0.05,
             p = c(1e-10, 1e-3, 1e-9, 0.5, 1e-12, 0.2)[seq_len(n)])
}

test_that("REGENIE binary-trait output is read, keeping the additive test", {
  r <- make_rows()
  regenie <- data.frame(CHROM = r$chr, GENPOS = r$pos, ID = r$id, ALLELE0 = "A", ALLELE1 = "G",
                        A1FREQ = r$af, INFO = 1, N = 1000, TEST = "ADD", BETA = r$beta,
                        SE = r$se, CHISQ = 10, LOG10P = -log10(r$p), EXTRA = NA)
  dom <- regenie
  dom$TEST <- "DOM"
  dom$LOG10P <- 50
  d <- standardise_sumstats(rbind(regenie, dom))
  expect_equal(nrow(d), nrow(regenie))
  expect_equal(d$logp, -log10(r$p))
  expect_equal(d$chr, c("1", "X", "1", "X", "1", "X"))
  expect_equal(attr(d, "columns")[["log10p"]], "LOG10P")
  expect_equal(attr(d, "columns")[["eaf"]], "A1FREQ")
})

test_that("SAIGE binary-trait output is read (p.value, not p.value.NA)", {
  r <- make_rows()
  saige <- data.frame(CHR = paste0("chr", sub("23", "X", r$chr)), POS = r$pos, SNPID = r$id,
                      Allele1 = "A", Allele2 = "G", AC_Allele2 = 10, AF_Allele2 = r$af,
                      imputationInfo = 1, N = 1000, BETA = r$beta, SE = r$se, Tstat = 1,
                      p.value = r$p, p.value.NA = 0.9, Is.SPA.converge = 1, varT = 1,
                      varTstar = 1, AF.Cases = 0.1, AF.Controls = 0.1, N.Cases = 100,
                      N.Controls = 900, check.names = FALSE)
  d <- standardise_sumstats(saige)
  expect_equal(attr(d, "columns")[["p"]], "p.value")
  expect_equal(attr(d, "columns")[["snp"]], "SNPID")
  expect_equal(d$logp, -log10(r$p))
  expect_equal(d$maf, pmin(r$af, 1 - r$af))
  expect_equal(d$beta, ifelse(r$af > 0.5, -r$beta, r$beta))
})

test_that("FinnGen-style output prefers the -log10 P column", {
  r <- make_rows()
  fg <- data.frame(`#chrom` = r$chr, pos = r$pos, ref = "A", alt = "G", rsids = r$id,
                   nearest_genes = "GENE", pval = 0, mlogp = -log10(r$p), beta = r$beta,
                   sebeta = r$se, af_alt = r$af, check.names = FALSE)
  d <- standardise_sumstats(fg)
  expect_equal(attr(d, "columns")[["log10p"]], "mlogp")
  expect_equal(d$logp, -log10(r$p))
})

test_that("a column of log10(P) rather than -log10(P) is turned round", {
  x <- data.frame(CHR = 1, BP = 1:3, LOG10P = log10(c(1e-8, 1e-3, 0.5)))
  msgs <- testthat::capture_messages(d <- standardise_sumstats(x, quiet = FALSE))
  expect_true(any(grepl("holds log10", msgs)))
  expect_equal(d$logp, c(8, 3, -log10(0.5)))
})

test_that("chromosome X is laid out after the autosomes and labelled", {
  gw <- simulate_gwas(n_variants = 2e4, seed = 2)
  x <- gw[gw$CHR <= 3, ]
  x$CHR[x$CHR == 3] <- 23
  p <- brownstone(x, effects = "none")
  sc <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x
  expect_equal(sc$get_labels(), c("1", "2", "X"))
})

test_that("a plot whose hits are all suggestive draws without warnings", {
  skip_if_not_installed("ragg")
  gw <- simulate_gwas(n_variants = 2e4, seed = 3)
  gw <- gw[gw$P > 5e-8, ]
  f <- tempfile(fileext = ".png")
  on.exit(unlink(f))
  expect_no_warning(ggplot2::ggsave(f, brownstone(gw, suggestive = 1e-4), width = 6, height = 4,
                                    device = ragg::agg_png))
})
