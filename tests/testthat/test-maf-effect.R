gw <- simulate_gwas(n_variants = 5e4, seed = 1)

point_layer <- function(plot, which = c("solid", "hollow", "cloud")) {
  which <- match.arg(which)
  i <- which(vapply(plot$layers, function(l) {
    if (!inherits(l$geom, "GeomPoint")) return(FALSE)
    switch(which,
           solid = identical(l$aes_params$colour, "white"),
           hollow = identical(l$aes_params$fill, "white"),
           cloud = identical(l$aes_params$colour, .ink$noise))
  }, logical(1)))
  ggplot2::layer_data(plot, i)
}
y_inverse <- function(p, y) {
  ggplot2::ggplot_build(p)$layout$panel_scales_y[[1]]$get_transformation()$inverse(y)
}

test_that("the effect-vs-frequency plot builds and renders on the default pdf() device", {
  p <- effect_vs_maf(gw)
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit({ grDevices::dev.off(); unlink(f) }, add = TRUE)
  expect_no_error(print(p))
})

test_that("hits sit at log10(MAF) and |beta|", {
  p <- effect_vs_maf(gw, labels = NULL)
  d <- standardise_sumstats(gw)
  hits <- d[d$logp >= -log10(5e-8), ]
  solid <- point_layer(p, "solid")
  expect_equal(nrow(solid), nrow(hits))
  expect_equal(sort(solid$x), sort(log10(hits$maf)))
  expect_equal(sort(y_inverse(p, solid$y)), sort(abs(hits$beta)))
  # Hits get one colour: frequency is already the x axis.
  expect_length(unique(solid$fill), 1)
  # The grey cloud is thinned.
  expect_lt(nrow(point_layer(p, "cloud")), sum(d$logp < -log10(5e-8)))
})

test_that("signed effects keep their sign, oriented to the minor allele", {
  p <- effect_vs_maf(gw, effects = "signed", labels = NULL)
  y <- y_inverse(p, point_layer(p, "solid")$y)
  expect_true(any(y < 0) && any(y > 0))
  expect_equal(p$labels$y, quote(beta ~ "(minor allele)"))
})

test_that("only genome-wide significant loci are labelled, up to max_labels", {
  repel <- function(p) {
    i <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomTextRepel"), logical(1)))
    lab <- ggplot2::layer_data(p, i)$label
    lab[lab != ""]
  }
  loci <- find_loci(gw)
  lab <- repel(effect_vs_maf(gw, suggestive = 1e-5))
  expect_length(lab, min(nrow(loci), 10))
  expect_equal(lab[1], loci$gene[1])
  expect_length(repel(effect_vs_maf(gw, max_labels = 2)), 2)
  expect_gt(nrow(point_layer(effect_vs_maf(gw, suggestive = 1e-5), "hollow")), 0)
})

test_that("frequency and effect columns are both required", {
  expect_error(effect_vs_maf(gw[, names(gw) != "EAF"]), "allele frequency")
  expect_error(suppressWarnings(effect_vs_maf(gw[, c("CHR", "BP", "P", "EAF")])), "effect size")
})

test_that("odds ratios are labelled as the odds ratio of the risk allele", {
  or <- gw
  or$OR <- exp(or$BETA)
  or$BETA <- NULL
  expect_equal(effect_vs_maf(or)$labels$y, "OR (risk allele)")
})

# A synthetic GWAS whose standard errors are exactly k / sqrt(2f(1 - f)): the
# detection boundary's fit is exact, and its curve stays inside the axis (the
# lead effects top the boundary's rare end).
boundary_gwas <- function(k = 0.02) {
  f <- c(0.4, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01)
  data.frame(CHR = "1", BP = seq_along(f) * 1000,
             P = c(1e-9, 1e-9, 1e-9, 1e-3, 0.5, 0.5, 0.5),
             BETA = c(2.0, -1.5, 1.0, 0, 0, 0, 0),
             SE = k / sqrt(2 * f * (1 - f)), EAF = f)
}

line_layer <- function(plot) {
  i <- which(vapply(plot$layers, function(l) inherits(l$geom, "GeomLine"), logical(1)))
  if (!length(i)) return(NULL)
  ggplot2::layer_data(plot, i[1])
}

test_that("the detection boundary is the study's own typical SE at each frequency", {
  d <- boundary_gwas()
  expect_equal(.se_unit(standardise_sumstats(d)), 0.02)
  p <- effect_vs_maf(d, labels = NULL)
  ld <- line_layer(p)
  z <- stats::qnorm(1 - 5e-8 / 2) + stats::qnorm(0.8)
  f <- 10^ld$x
  expect_equal(y_inverse(p, ld$y), z * 0.02 / sqrt(2 * f * (1 - f)))
})

test_that("the boundary tracks the median SE of variants at each frequency", {
  p <- effect_vs_maf(gw, labels = NULL)
  ld <- line_layer(p)
  d <- standardise_sumstats(gw)
  z <- stats::qnorm(1 - 5e-8 / 2) + stats::qnorm(0.8)
  mid <- d$maf >= 0.1 & d$maf < 0.2
  expect_gt(sum(mid), 1000)
  expect_equal(approx(ld$x, y_inverse(p, ld$y), xout = log10(stats::median(d$maf[mid])))$y,
               z * stats::median(d$se[mid]), tolerance = 0.1)
})

test_that("the boundary is labelled with its power at the common end", {
  p <- effect_vs_maf(gw, labels = NULL)
  i <- which(vapply(p$layers, function(l) {
    inherits(l$geom, "GeomText") && grepl("power", paste(l$aes_params$label, collapse = ""))
  }, logical(1)))
  expect_length(i, 1)
  ld <- ggplot2::layer_data(p, i)
  expect_equal(ld$label, "80% power")
  expect_equal(ld$x, log10(0.5))
  z <- stats::qnorm(1 - 5e-8 / 2) + stats::qnorm(0.8)
  expect_equal(y_inverse(p, ld$y),
               z * .se_unit(standardise_sumstats(gw)) / sqrt(2 * 0.5 * 0.5))
})

test_that("signed effects mirror the boundary; several powers draw several curves", {
  ps <- effect_vs_maf(boundary_gwas(), effects = "signed", labels = NULL)
  ld <- line_layer(ps)
  expect_equal(unique(ld$group), c(1, 2))
  expect_equal(y_inverse(ps, ld$y[ld$group == 1]), -y_inverse(ps, ld$y[ld$group == 2]))
  pp <- effect_vs_maf(boundary_gwas(), power = c(0.5, 0.8), labels = NULL)
  l2 <- line_layer(pp)
  expect_equal(unique(l2$group), c(1, 2))
  zq <- stats::qnorm(1 - 5e-8 / 2)
  r <- y_inverse(pp, l2$y[l2$group == 2]) / y_inverse(pp, l2$y[l2$group == 1])
  expect_equal(r, rep((zq + stats::qnorm(0.8)) / (zq + stats::qnorm(0.5)), 200))
})

test_that("the boundary is optional, and absent when no SE is estimable", {
  p <- effect_vs_maf(gw, detection = FALSE, labels = NULL)
  expect_null(line_layer(p))
  expect_null(.detection_boundary(NA_real_, 1e-4, 5e-8, 0.8, FALSE))
  expect_null(.boundary_layers(NULL, 10))
})

test_that("power is a level, not a percentage", {
  expect_error(effect_vs_maf(gw, power = 80), "between 0 and 1")
})
