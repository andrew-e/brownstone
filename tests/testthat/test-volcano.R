gw <- simulate_gwas(n_variants = 5e4, seed = 1)

# Point layers are told apart by their fixed parameters: solid (genome-wide
# significant) points have a white outline, hollow (suggestive) ones a white
# fill, and the background cloud neither.
point_layer <- function(plot, which = c("solid", "hollow", "cloud")) {
  which <- match.arg(which)
  i <- which(vapply(plot$layers, function(l) {
    if (!inherits(l$geom, "GeomPoint")) return(FALSE)
    switch(which,
           solid = identical(l$aes_params$colour, "white"),
           hollow = identical(l$aes_params$fill, "white"),
           cloud = is.null(l$aes_params$fill) && !identical(l$aes_params$colour, "white"))
  }, logical(1)))
  ggplot2::layer_data(plot, i)
}
x_inverse <- function(p, x) {
  ggplot2::ggplot_build(p)$layout$panel_scales_x[[1]]$get_transformation()$inverse(x)
}
repel_labels <- function(p) {
  i <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomTextRepel"), logical(1)))
  lab <- ggplot2::layer_data(p, i)$label
  lab[lab != ""]
}

test_that("the volcano plot builds and renders on the default pdf() device", {
  p <- volcano_plot(gw)
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit({ grDevices::dev.off(); unlink(f) }, add = TRUE)
  expect_no_error(print(p))
  expect_no_error(print(volcano_plot(gw, suggestive = 1e-5)))
})

test_that("hits sit at their signed effect and -log10(P)", {
  p <- volcano_plot(gw, labels = NULL)
  d <- standardise_sumstats(gw)
  hits <- d[d$logp >= -log10(5e-8), ]
  solid <- point_layer(p, "solid")
  expect_equal(nrow(solid), nrow(hits))
  expect_equal(sort(x_inverse(p, solid$x)), sort(hits$beta))
  expect_true(any(solid$x < 0) && any(solid$x > 0))
  below_break <- hits$logp < 20
  expect_equal(sort(solid$y[solid$y < 20]), sort(hits$logp[below_break]))
})

test_that("the null cloud is grey and thinned", {
  p <- volcano_plot(gw, labels = NULL)
  cloud <- point_layer(p, "cloud")
  expect_gt(nrow(cloud), 0)
  expect_lt(nrow(cloud), sum(gw$P >= 5e-8))
  expect_equal(unique(cloud$colour), .ink$noise)
})

test_that("loci are labelled once each, most significant first", {
  lab <- repel_labels(volcano_plot(gw))
  loci <- find_loci(gw)
  expect_length(lab, min(nrow(loci), 20))
  expect_equal(lab[1], loci$gene[1])
  expect_length(repel_labels(volcano_plot(gw, max_labels = 3)), 3)
  expect_length(repel_labels(volcano_plot(gw, labels = "snp", max_labels = 2)), 2)
})

test_that("suggestive peaks are rings; the suggestive line must be weaker", {
  p <- volcano_plot(gw, suggestive = 1e-5, labels = NULL)
  expect_gt(nrow(point_layer(p, "hollow")), 0)
  expect_error(volcano_plot(gw, suggestive = 1e-9), "larger P-value")
})

test_that("odds ratios are drawn on the log scale and labelled as odds ratios", {
  or <- gw
  or$OR <- exp(or$BETA)
  or$BETA <- NULL
  p <- volcano_plot(or, labels = NULL)
  expect_equal(p$labels$x, "OR (minor allele)")
  d <- standardise_sumstats(or)   # log(OR), oriented to the minor allele
  expect_equal(sort(x_inverse(p, point_layer(p, "solid")$x)),
               sort(d$beta[d$logp >= -log10(5e-8)]))
})

test_that("the effect axis breaks only for outlying estimates", {
  scale_name <- function(p) ggplot2::ggplot_build(p)$layout$panel_scales_x[[1]]$get_transformation()$name
  out <- gw
  out$BETA[which(out$P > 0.5)[1]] <- 50
  expect_equal(scale_name(volcano_plot(out, labels = NULL)), "broken_log")
  expect_equal(scale_name(volcano_plot(out, labels = NULL, effect_scale = "linear")), "identity")
})

test_that("without effects there is no volcano; without frequencies, one colour", {
  expect_error(suppressWarnings(volcano_plot(gw[, c("CHR", "BP", "P")])), "effect size")
  expect_message(p <- volcano_plot(gw[, names(gw) != "EAF"], labels = NULL), "one colour")
  expect_length(unique(point_layer(p, "solid")$fill), 1)
  expect_equal(unique(point_layer(p, "cloud")$colour), .ink$noise)
})
