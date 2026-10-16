toy <- function() {
  data.frame(
    CHR = c(1, 1, 1, 1, 2, 2),
    BP = c(1e6, 1.2e6, 1.9e6, 5e6, 3e6, 3.1e6),
    P = c(1e-10, 1e-9, 1e-12, 1e-3, 1e-20, 1e-9)
  )
}

test_that("distance clumping finds the expected loci", {
  l <- find_loci(toy(), annotate = FALSE)
  expect_equal(nrow(l), 3)
  expect_equal(l$chr, c("2", "1", "1"))
  expect_equal(l$pos, c(3e6, 1.9e6, 1e6))
  expect_equal(l$n_sig, c(2, 1, 2))
  expect_equal(l$start[3], 1e6)
  expect_equal(l$end[3], 1.2e6)
  # No ID column: lead variants are named chr:pos.
  expect_equal(l$snp, c("2:3000000", "1:1900000", "1:1000000"))
})

test_that("a wider window merges nearby loci", {
  l <- find_loci(toy(), window = 1e6, annotate = FALSE)
  expect_equal(nrow(l), 2)
})

test_that("no significant variants gives an empty table", {
  l <- find_loci(data.frame(CHR = 1, BP = 1, P = 0.5), annotate = FALSE)
  expect_equal(nrow(l), 0)
})

test_that("loci are annotated with the nearest gene", {
  # PCSK9 (chr1) and APOE (chr19) in both builds.
  l38 <- find_loci(data.frame(CHR = c(1, 19), BP = c(55043980, 44905000), P = 1e-10))
  expect_setequal(l38$gene, c("PCSK9", "APOE"))
  l37 <- find_loci(data.frame(CHR = c(1, 19), BP = c(55510000, 45411000), P = 1e-10),
                   build = "GRCh37")
  expect_setequal(l37$gene, c("PCSK9", "APOE"))
})

test_that("nearest_gene handles inside, between and custom genes", {
  g <- data.frame(chr = "1", start = c(100, 300), end = c(200, 400),
                  strand = c("+", "-"), gene = c("A", "B"))
  ng <- nearest_gene(c(1, 1, 1, 2), c(150, 260, 1000, 50), genes = g)
  expect_equal(ng$gene, c("A", "B", "B", NA))
  expect_equal(ng$distance, c(0, 40, 600, NA))
  expect_equal(nearest_gene(1, 1000, genes = g, max_distance = 100)$gene, NA_character_)
})
