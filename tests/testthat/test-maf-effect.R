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
