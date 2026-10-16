test_that("columns from common GWAS tools are detected", {
  formats <- list(
    plink2 = data.frame(`#CHROM` = 1, POS = 100, ID = "rs1", A1_FREQ = 0.2,
                        BETA = 0.1, SE = 0.02, P = 1e-6, check.names = FALSE),
    regenie = data.frame(CHROM = 1, GENPOS = 100, ID = "rs1", A1FREQ = 0.2,
                         BETA = 0.1, SE = 0.02, LOG10P = 6),
    bolt = data.frame(SNP = "rs1", CHR = 1, BP = 100, A1FREQ = 0.2, BETA = 0.1,
                      SE = 0.02, P_BOLT_LMM = 1e-6),
    gwas_ssf = data.frame(chromosome = 1, base_pair_location = 100,
                          effect_allele_frequency = 0.2, beta = 0.1,
                          standard_error = 0.02, p_value = 1e-6, rsid = "rs1")
  )
  for (nm in names(formats)) {
    d <- standardise_sumstats(formats[[nm]])
    expect_equal(d$chr, "1", info = nm)
    expect_equal(d$pos, 100, info = nm)
    expect_equal(d$snp, "rs1", info = nm)
    expect_equal(d$logp, 6, tolerance = 1e-9, info = nm)
    expect_equal(d$beta, 0.1, info = nm)
    expect_equal(d$se, 0.02, info = nm)
    expect_equal(d$maf, 0.2, info = nm)
  }
})

test_that("explicit column names override detection", {
  x <- data.frame(c = "chr2", where = 5, pv = 0.01)
  expect_error(standardise_sumstats(x), "chromosome and position")
  d <- standardise_sumstats(x, chr = "c", pos = "where", p = "pv")
  expect_equal(d$chr, "2")
  expect_equal(d$logp, 2)
  expect_error(standardise_sumstats(x, chr = "nope", pos = "where", p = "pv"), "not in the data")
})

test_that("chromosome names are normalised", {
  x <- data.frame(CHR = c("chr1", "23", "chrX", "24", "26", "M", "6_hap"),
                  BP = 1:7, P = 0.5)
  d <- suppressMessages(standardise_sumstats(x))
  expect_equal(d$chr, c("1", "X", "X", "Y", "MT", "MT"))
})

test_that("frequencies are folded and effects oriented to the minor allele", {
  x <- data.frame(CHR = 1, BP = 1:2, P = 0.1, BETA = c(0.3, 0.3), SE = 0.1,
                  EAF = c(0.1, 0.9))
  d <- standardise_sumstats(x)
  expect_equal(d$maf, c(0.1, 0.1))
  expect_equal(d$beta, c(0.3, -0.3))
  expect_equal(attr(d, "effect_allele"), "minor")
  d2 <- standardise_sumstats(x, effect_allele = "as_is")
  expect_equal(d2$beta, c(0.3, 0.3))
  expect_equal(attr(d2, "effect_allele"), "as_is")
})

test_that("odds ratios are used on the log scale", {
  x <- data.frame(CHR = 1, BP = 1, P = 0.1, OR = 2, SE = 0.1)
  d <- standardise_sumstats(x)
  expect_equal(d$beta, log(2))
  expect_equal(attr(d, "effect_scale"), "log_or")
})

test_that("P-values that underflowed to zero are recovered from beta/se", {
  x <- data.frame(CHR = 1, BP = 1:2, P = c(0, 0), BETA = c(4, 0.1), SE = c(0.1, NA))
  d <- suppressWarnings(standardise_sumstats(x))
  expect_gt(d$logp[1], 340)   # z = 40
  expect_lt(d$logp[1], 352)
  expect_warning(standardise_sumstats(x), "could not be recovered")
})

test_that("P-values written as text are read exactly below double precision", {
  x <- data.frame(CHR = 1, BP = 1:3, P = c("1.24e-652", "0.5", "3E-10"),
                  BETA = 0.1, SE = 0.01)
  d <- standardise_sumstats(x)
  expect_equal(d$logp, c(652 - log10(1.24), -log10(0.5), 10 - log10(3)))
})

test_that("bad P-values are caught", {
  x <- data.frame(CHR = 1, BP = 1, P = 7)
  expect_error(standardise_sumstats(x), "outside \\[0, 1\\]")
})

test_that("extra columns can be kept", {
  x <- data.frame(CHR = 1, BP = 1, P = 0.1, GENE = "ABC")
  expect_equal(standardise_sumstats(x, keep = "GENE")$GENE, "ABC")
  expect_error(standardise_sumstats(x, keep = "NOPE"), "not in the data")
})

test_that("variants without an ID column are named chr:pos", {
  x <- data.frame(CHR = c(1, 23), BP = c(55000000, 1234), P = 0.1)
  expect_equal(standardise_sumstats(x)$snp, c("1:55000000", "X:1234"))
  # Internally the IDs are only built for the rows that need them.
  d <- .standardise(x, ids = FALSE)
  expect_true(all(is.na(d$snp)))
  expect_equal(.variant_id(d, 2:1), c("X:1234", "1:55000000"))
  # An ID column is used as it is, missing values included.
  y <- data.frame(CHR = 1, BP = 1:2, P = 0.1, SNP = c("rs1", NA))
  expect_equal(.variant_id(.standardise(y, ids = FALSE), 1:2), c("rs1", NA))
})
