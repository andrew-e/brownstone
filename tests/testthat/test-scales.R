test_that("y compression is linear below the break and monotone above", {
  comp <- .y_compressor(300, sig = 7.3, y_compress = 60)
  expect_equal(comp$at, 60)
  y <- c(0, 10, 59.9)
  expect_equal(comp$f(y), y)
  up <- comp$f(c(61, 100, 200, 300))
  expect_true(all(diff(up) > 0))
})

test_that("the compressed zone gets a fixed share per decade, within limits", {
  share <- function(max_logp) {
    comp <- .y_compressor(max_logp, sig = 7.3)
    1 - comp$at / comp$f(max_logp)
  }
  expect_equal(share(200), 0.15)                  # one decade above 20
  expect_equal(share(41), 0.12)                   # floor (just past the 2x margin)
  expect_equal(share(1e5), 0.35)                  # ceiling
})

test_that("the automatic break sits at 2.5 times the significance line", {
  expect_equal(.auto_compress_at(-log10(5e-8)), 20)
  expect_equal(.auto_compress_at(-log10(1e-5)), 15)
  expect_null(.y_compressor(19.9, sig = 7.3)$at)    # nothing above the break
  expect_null(.y_compressor(21, sig = 7.3)$at)      # too little above it (needs ~2x)
  expect_null(.y_compressor(25, sig = 7.3)$at)
  expect_null(.y_compressor(30, sig = 7.3)$at)      # still under 2 x 20
  expect_equal(.y_compressor(41, sig = 7.3)$at, 20) # past 2 x 20, the break kicks in
  expect_equal(.y_compressor(21, sig = 7.3, y_compress = 20)$at, 20)   # asked for
  expect_null(.y_compressor(300, sig = 7.3, FALSE)$at)
})

test_that("the compressed zone always has a label above the break", {
  comp <- .y_compressor(50, sig = 7.3)
  expect_equal(comp$at, 20)
  expect_equal(max(.y_breaks(comp, 50)$labels), 50)
  # A tail long enough for 1-2-5 values gets those, as before.
  expect_equal(.y_breaks(.y_compressor(300, sig = 7.3), 300)$labels[-(1:5)], c(50, 100, 200))
})

test_that("a tail just above the break is never stretched", {
  comp <- .y_compressor(21, sig = 7.3, y_compress = 20)
  # share of the axis given to 20-21 is at most its linear share (1/21)
  expect_lte(1 - comp$at / comp$f(21), 1 / 21 + 1e-9)
})

test_that("the effect axis breaks only for outlying estimates", {
  bulk <- c(0.02, 0.03, 0.04, 0.05, 0.06, 0.08)
  expect_equal(.effect_axis(bulk, top = 0.1)$scale, "linear")
  ax <- .effect_axis(c(bulk, 1.3), top = 1.6)
  expect_equal(ax$scale, "broken")
  expect_equal(ax$at, 0.2)                          # 3 x median 0.055, rounded up
  expect_equal(ax$trans$transform(0.1), 0.1)        # linear below the break
  expect_equal(ax$trans$inverse(ax$trans$transform(c(0.5, 1.3))), c(0.5, 1.3))
  expect_equal(ax$trans$transform(-0.5), -ax$trans$transform(0.5))
  # wide intervals alone don't trigger a break
  expect_equal(.effect_axis(bulk, top = 1.6, max_estimate = 0.08)$scale, "linear")
  expect_equal(.effect_axis(bulk, top = 1, scale = "sqrt")$scale, "sqrt")
  b <- .track_breaks(ax, top = 1.6, signed = FALSE, log_or = FALSE)
  expect_true(all(c(0, 0.2, 1) %in% b))
})

test_that("the compressed zone of the effect axis always has a tick", {
  bulk <- c(0.02, 0.03, 0.04, 0.05, 0.06, 0.08)
  # 0.2 to 0.62: no 1-2-5 value half a decade above the break fits.
  ax <- .effect_axis(bulk, top = 0.62, max_estimate = 0.62)
  expect_equal(ax$at, 0.2)
  b <- .track_breaks(ax, top = 0.62, signed = FALSE, log_or = FALSE)
  expect_equal(b, c(0, 0.1, 0.2, 0.5))
  expect_equal(.track_breaks(ax, top = 0.62, signed = TRUE, log_or = FALSE),
               c(-0.5, -0.1, 0, 0.1, 0.5))
  # A full-size panel keeps the tick at the break.
  expect_equal(.track_breaks(ax, top = 0.62, signed = TRUE, log_or = FALSE, roomy = TRUE),
               c(-0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5))
})

test_that("ticks above the break are spaced on the log scale", {
  comp <- .y_compressor(317, sig = 7.3)
  expect_equal(.y_breaks(comp, 317)$labels, c(0, 5, 10, 15, 20, 50, 100, 200))
})

test_that("effect breaks are symmetric round numbers including zero", {
  expect_equal(.effect_breaks(c(-0.45, 0.45), "sqrt"), c(-0.4, -0.1, 0, 0.1, 0.4))
  expect_equal(.effect_breaks(c(-0.6, 0.6), "sqrt"), c(-0.5, -0.1, 0, 0.1, 0.5))
  expect_equal(.effect_breaks(c(-0.3, 0.3), "linear"), c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3))
  b <- .effect_breaks(log(c(1 / 1.6, 1.6)), "sqrt", log_or = TRUE)
  expect_equal(b, log(c(1 / 1.5, 1 / 1.1, 1, 1.1, 1.5)))
  expect_equal(.effect_labels(c(-0.1, 0, 0.5)), c("-0.1", "0", "0.5"))
  expect_equal(.effect_labels(log(c(1, 2)), log_or = TRUE), c("1", "2"))
})

test_that("MAF classes are ordered from common to rare", {
  f <- .maf_classes(c(0.3, 0.02, 0.001, 0.05, 0.01), c(0.01, 0.05))
  expect_equal(levels(f), c(">= 5%", "1-5%", "< 1%"))
  expect_equal(as.character(f), c(">= 5%", "1-5%", "< 1%", ">= 5%", "1-5%"))
  lab <- .maf_legend_labels(levels(f))
  expect_true(is.expression(lab))
  expect_identical(lab[[1]], quote("" >= "5%"))
})

test_that("the LD ceiling from allele frequencies is right", {
  expect_equal(.r_max(0.3, 0.3), 1)
  expect_equal(.r_max(0.01, 0.3), .r_max(0.3, 0.01))
  expect_equal(.r_max(0.1, 0.5), sqrt(0.1 * 0.5 / (0.5 * 0.9)))
})

test_that("the P-value plotmath reads naturally", {
  expect_equal(.p_plotmath(5e-8), "5 %*% 10^-8")
  expect_equal(.p_plotmath(1e-5), "10^-5")
})
