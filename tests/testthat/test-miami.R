gw <- simulate_gwas(n_variants = 5e4, seed = 1)
gw2 <- replicate_gwas(gw)

# Every point layer of a panel, stacked.
all_points <- function(p) {
  i <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1)))
  do.call(rbind, lapply(i, function(k) {
    l <- ggplot2::layer_data(p, k)
    if (nrow(l)) l[, c("x", "y")]
  }))
}
solid_points <- function(p) {
  i <- which(vapply(p$layers, function(l) {
    inherits(l$geom, "GeomPoint") && identical(l$aes_params$colour, "white")
  }, logical(1)))
  ggplot2::layer_data(p, i)
}
y_range <- function(p) ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y.range
x_range <- function(p) ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x.range
# Each half is named in the margin: above (a panel title, or for the upper
# -log10(P) panel the top x-axis title, below the legend), or below (a caption).
tags <- function(p) {
  top <- p$scales$get_scales("x")$secondary.axis$name
  c(unlist(p$labels[c("title", "caption")]), if (is.character(top)) c(x.top = top))
}

test_that("the Miami plot builds and renders on the default pdf() device", {
  p <- miami_plot(gw, gw2)
  expect_s3_class(p, "patchwork")
  expect_length(p$patches$plots, 3)   # plus the last panel: four in all
  for (k in 1:4) expect_silent(ggplot2::ggplot_build(p[[k]]))
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit({ grDevices::dev.off(); unlink(f) }, add = TRUE)
  expect_no_error(print(p))
  expect_no_error(print(miami_plot(gw, gw2, suggestive = 1e-5, label_angle = 90,
                                   maf_scale = "continuous")))
})

test_that("the first GWAS rises, the second hangs, on one -log10(P) scale", {
  p <- miami_plot(gw, gw2)
  expect_true(all(all_points(p[[1]])$y >= 0))
  expect_true(all(all_points(p[[2]])$y <= 0))
  expect_equal(y_range(p[[1]]), -rev(y_range(p[[2]])))
  # The scale is set by the stronger GWAS, whichever half it is in.
  weak <- replicate_gwas(gw, shrink = 0.3)
  strong_top <- miami_plot(gw, weak, effects = FALSE)
  strong_bottom <- miami_plot(weak, gw, effects = FALSE)
  expect_equal(y_range(strong_top[[1]]), y_range(strong_bottom[[1]]))
  expect_equal(y_range(strong_top[[2]]), y_range(strong_bottom[[2]]))
})

test_that("one bold rule separates the -log10(P) plot from the |beta| plot", {
  p <- miami_plot(gw, gw2)
  is_bold <- function(l) inherits(l$geom, "GeomHline") && identical(l$aes_params$linewidth, .ground_width)
  expect_true(any(vapply(p[[3]]$layers, is_bold, logical(1))))
  # Nothing else is drawn bold: not the axis between the -log10(P) halves...
  for (k in c(1, 2, 4)) {
    expect_false(inherits(p[[k]]$theme$axis.line.x.bottom, "element_line"))
    expect_false(any(vapply(p[[k]]$layers, is_bold, logical(1))))
  }
  # ...and the halves share one -log10(P) axis title.
  expect_equal(p$patches$layout$axis_titles, "collect_y")
})

test_that("both halves share one genome axis, even when one lacks a chromosome", {
  no22 <- gw2[gw2$CHR != 22, ]
  p <- miami_plot(gw, no22)
  expect_equal(x_range(p[[1]]), x_range(p[[2]]))
  expect_equal(x_range(p[[1]]), x_range(p[[4]]))
  expect_equal(x_range(p[[1]]), x_range(miami_plot(gw, gw2)[[1]]))
})

test_that("the effect track shows |beta| of the first GWAS above zero, the second below", {
  p <- miami_plot(gw, gw2)
  # Layer data are on the axis's (possibly broken) scale.
  y_scale <- ggplot2::ggplot_build(p[[4]])$layout$panel_scales_y[[1]]
  y <- y_scale$get_transformation()$inverse(solid_points(p[[4]])$y)
  up <- abs(find_loci(gw, annotate = FALSE)$beta)
  down <- abs(find_loci(gw2, annotate = FALSE)$beta)
  expect_true(all(round(up, 8) %in% round(y[y > 0], 8)))
  expect_true(all(round(down, 8) %in% round(-y[y < 0], 8)))
  # Tick labels are magnitudes on both sides.
  labs <- ggplot2::ggplot_build(p[[4]])$layout$panel_params[[1]]$y$get_labels()
  expect_false(any(grepl("-", vapply(labs, deparse, ""), fixed = TRUE)))
})

test_that("one palette for both halves, with an Unknown class from either", {
  unk <- gw2
  unk$EAF[unk$P < 1e-6 & unk$CHR == 19] <- NA
  p <- miami_plot(gw, unk)
  fill <- ggplot2::ggplot_build(p[[1]])$plot$scales$get_scales("fill")
  expect_true("Unknown" %in% fill$get_limits())
})

test_that("without comparable effects there is no track, with a message", {
  or <- gw2
  or$OR <- exp(or$BETA)
  or$BETA <- NULL
  expect_message(p <- miami_plot(gw, or), "odds ratios")
  expect_length(p$patches$plots, 1)   # two panels
  expect_message(p <- suppressWarnings(miami_plot(gw, gw2[, c("CHR", "BP", "P")])),
                 "effect track")
  expect_length(p$patches$plots, 1)
  expect_length(miami_plot(gw, gw2, effects = FALSE)$patches$plots, 1)
})

test_that("halves are named, and labels can come from a table", {
  p <- miami_plot(gw, gw2)
  expect_equal(tags(p[[1]]), c(x.top = "gw"))
  expect_equal(tags(p[[2]]), c(caption = "gw2"))
  expect_equal(tags(p[[4]]), c(title = "gw", caption = "gw2"))
  p <- miami_plot(gw, gw2, names = c("Men", "Women"))
  expect_equal(tags(p[[2]]), c(caption = "Women"))
  expect_error(miami_plot(gw, gw2, names = "Men"), "names")

  # A label for a variant in only one GWAS is drawn there, without a warning.
  label_text <- function(panel) {
    i <- which(vapply(panel$layers, function(l) inherits(l$geom, "GeomTextRepel"), logical(1)))
    lab <- ggplot2::layer_data(panel, i)$label
    lab[lab != ""]
  }
  k <- which(gw$P > 0.5)[1]
  only_top <- gw2[gw2$SNP != gw$SNP[k], ]
  lab <- data.frame(snp = gw$SNP[k], label = "QUIET")
  expect_no_warning(p <- miami_plot(gw, only_top, labels = lab))
  expect_equal(label_text(p[[1]]), "QUIET")
  expect_warning(miami_plot(gw, gw2, labels = data.frame(snp = "nowhere", label = "X")),
                 "either GWAS")
})
