gw <- simulate_gwas(n_variants = 5e4, seed = 1)

# Point layers are told apart by their fixed parameters: solid (genome-wide
# significant) points have a white outline, hollow (suggestive) ones a white
# fill.
point_layer <- function(plot, which = c("solid", "hollow")) {
  which <- match.arg(which)
  param <- if (which == "solid") "colour" else "fill"
  i <- which(vapply(plot$layers, function(l) {
    inherits(l$geom, "GeomPoint") && identical(l$aes_params[[param]], "white")
  }, logical(1)))
  ggplot2::layer_data(plot, i)
}

test_that("the full plot builds, with an effect track", {
  p <- manhattan(gw)
  expect_s3_class(p, "patchwork")
  expect_silent(ggplot2::ggplot_build(p[[1]]))
  expect_silent(ggplot2::ggplot_build(p[[2]]))
})

test_that("the plot renders on the default pdf() device", {
  # Raw Unicode text (beta, minus, >=) can't be drawn by pdf() in a single-byte
  # encoding; everything special must go through plotmath.
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit({ grDevices::dev.off(); unlink(f) }, add = TRUE)
  or <- gw
  or$OR <- exp(or$BETA)
  or$BETA <- NULL
  expect_no_error(print(manhattan(gw, suggestive = 1e-5)))
  expect_no_error(print(manhattan(or)))
})

test_that("a broken axis is a white gap, with no rules or marks along it", {
  p <- manhattan(gw)
  for (k in 1:2) {
    geoms <- vapply(p[[k]]$layers, function(l) class(l$geom)[1], "")
    expect_false("GeomCustomAnn" %in% geoms)
  }
  # No rules along the gap: the main panel's only rule is the significance line.
  geoms <- vapply(p[[1]]$layers, function(l) class(l$geom)[1], "")
  expect_false("GeomSegment" %in% geoms)
  expect_equal(sum(geoms == "GeomHline"), 1)
  mark <- .break_mark(1, 2, "y")
  expect_s3_class(mark$geom, "GeomRect")
  expect_equal(mark$aes_params$fill, "white")
  # Signed effects break on both sides of zero.
  ax <- list(at = 0.2, trans = scales::transform_identity())
  expect_length(.effect_break(ax, 1, signed = TRUE, axis = "x"), 2)
  expect_null(.effect_break(list(at = NULL), 1, signed = TRUE))
})

test_that("without effect columns the plot is a single ggplot", {
  # Without beta/se the underflowed APOE P-value can't be recovered: that
  # warning is expected here.
  suppressWarnings(expect_message(p <- manhattan(gw[, c("CHR", "BP", "P")]), "effect track"))
  expect_s3_class(p, "ggplot")
  expect_false(inherits(p, "patchwork"))
  expect_s3_class(manhattan(gw, effects = FALSE), "ggplot")
})

test_that("labels can come from a column or a table", {
  gw$GENE <- paste0("G", seq_len(nrow(gw)))
  expect_s3_class(manhattan(gw, labels = "GENE"), "patchwork")
  lead <- find_loci(gw)
  lab <- data.frame(snp = lead$snp, label = toupper(lead$gene),
                    fontface = rep(c("bold", "plain"), length.out = nrow(lead)))
  expect_s3_class(manhattan(gw, labels = lab), "patchwork")
  expect_s3_class(manhattan(gw, labels = NULL), "patchwork")
  expect_error(manhattan(gw, labels = data.frame(snp = "x")), "label")
})

test_that("a table of labels can name any variant, by ID or by position", {
  label_text <- function(p) {
    i <- which(vapply(p[[1]]$layers, function(l) inherits(l$geom, "GeomTextRepel"), logical(1)))
    lab <- ggplot2::layer_data(p[[1]], i)$label
    lab[lab != ""]
  }
  # A variant far from significance, at a position that as.character() writes
  # as "5.5e+07".
  k <- which(gw$CHR == 1 & gw$P > 0.5)[1]
  gw$BP[k] <- 55000000
  by_id <- data.frame(snp = gw$SNP[k], label = "QUIET")
  expect_equal(label_text(manhattan(gw, labels = by_id)), "QUIET")
  by_pos <- data.frame(chr = 1L, pos = 55000000L, label = "QUIET")
  expect_equal(label_text(manhattan(gw, labels = by_pos)), "QUIET")

  # Without an ID column, variants are named chr:pos.
  noid <- gw[, names(gw) != "SNP"]
  by_chrpos <- data.frame(snp = "1:55000000", label = "QUIET")
  expect_equal(label_text(manhattan(noid, labels = by_chrpos)), "QUIET")
  snp_labels <- label_text(manhattan(noid, labels = "snp"))
  top <- find_loci(noid, annotate = FALSE)[1, ]
  expect_equal(snp_labels[1], paste0(top$chr, ":", format(top$pos, scientific = FALSE)))
  expect_true(all(grepl("^[0-9XY]+:[0-9]+$", snp_labels)))
})

test_that("the effect track shows magnitudes by default, signs on request", {
  track_y <- function(p) point_layer(p[[2]], "solid")$y
  mag <- manhattan(gw)
  expect_true(all(track_y(mag) >= 0))
  signed <- manhattan(gw, effects = "signed")
  expect_true(any(track_y(signed) < 0))
  expect_s3_class(manhattan(gw, effects = TRUE), "patchwork")
  expect_s3_class(manhattan(gw, effects = "none"), "ggplot")
  expect_s3_class(manhattan(gw, track_direction = "down"), "patchwork")
})

test_that("suggestive peaks are drawn, labelled and kept off the track", {
  p <- manhattan(gw, suggestive = 1e-5)
  main <- ggplot2::ggplot_build(p[[1]])
  shapes <- main$plot$scales$get_scales("shape")
  expect_equal(shapes$get_limits(), c("gws", "suggestive"))

  # Hollow rings only for suggestive-only loci, never the base of a GWS tower.
  d <- standardise_sumstats(gw)
  gws_lead <- find_loci(gw, annotate = FALSE)
  weak <- point_layer(p[[1]], "hollow")
  lay <- .genome_layout(d)
  weak_pos <- weak$x - lay$offset[findInterval(weak$x, lay$start)]
  weak_chr <- lay$chr[findInterval(weak$x, lay$start)]
  near_gws <- vapply(seq_along(weak_pos), function(i) {
    any(gws_lead$chr == weak_chr[i] & abs(gws_lead$pos - weak_pos[i]) <= 5e5)
  }, logical(1))
  expect_gt(nrow(weak), 0)
  expect_false(any(near_gws))

  # The track has the GWS leads (solid) and the suggestive leads (hollow)...
  n_track <- nrow(point_layer(p[[2]], "solid"))
  expect_gte(n_track, nrow(gws_lead))
  n_sugg <- nrow(find_loci(gw, threshold = 1e-5, annotate = FALSE)) - nrow(gws_lead)
  expect_equal(nrow(point_layer(p[[2]], "hollow")), n_sugg)
  # ...unless they are kept off it.
  off <- manhattan(gw, suggestive = 1e-5, track_suggestive = FALSE)
  expect_equal(nrow(point_layer(off[[2]], "hollow")), 0)

  # Suggestive labels may sit between the two lines, beside their peaks,
  # rather than all above the significance line; without them, labels stay
  # above it.
  ylim_low <- function(p) {
    Filter(function(l) inherits(l$geom, "GeomTextRepel"), p[[1]]$layers)[[1]]$geom_params$ylim[1]
  }
  expect_lt(ylim_low(p), -log10(5e-8))
  expect_gt(ylim_low(p), 5)
  expect_gt(ylim_low(manhattan(gw)), -log10(5e-8))

  expect_error(manhattan(gw, suggestive = 1e-9), "larger P-value")
})

test_that("the track can show every hit, not just the leads", {
  cloud_layers <- function(p) {
    Filter(function(l) inherits(l$geom, "GeomPoint") && is.null(l$aes_params$fill) &&
             !identical(l$aes_params$colour, "white"), p[[2]]$layers)
  }
  leads <- manhattan(gw, track_variants = "leads")
  expect_true(all(vapply(cloud_layers(leads), function(l) nrow(l$data), 0) == 0))
  # signed effects always show lead variants only
  signed <- manhattan(gw, effects = "signed")
  expect_true(all(vapply(cloud_layers(signed), function(l) nrow(l$data), 0) == 0))

  all_gws <- manhattan(gw)                 # "all" is the default
  n_cloud <- sum(vapply(seq_along(all_gws[[2]]$layers), function(i) {
    l <- all_gws[[2]]$layers[[i]]
    if (inherits(l$geom, "GeomPoint") && isTRUE(l$aes_params$shape %in% c(1, 16))) {
      nrow(ggplot2::layer_data(all_gws[[2]], i))
    } else 0
  }, 0))
  n_sig <- sum(gw$P < 5e-8 & is.finite(gw$BETA))
  expect_gt(n_cloud, 0)
  expect_lt(n_cloud, n_sig)                    # leads are drawn separately
  expect_true(all(point_layer(all_gws[[2]], "solid")$y >= 0))

  with_sugg <- manhattan(gw, track_variants = "all", suggestive = 1e-5)
  rings <- Filter(function(l) identical(l$aes_params$shape, 1), with_sugg[[2]]$layers)
  expect_gt(nrow(rings[[1]]$data), 0)
})

test_that("there are four frequency classes by default", {
  p <- manhattan(gw)
  fill <- ggplot2::ggplot_build(p[[1]])$plot$scales$get_scales("fill")
  expect_equal(fill$get_limits(), c(">= 5%", "1-5%", "0.1-1%", "< 0.1%"))
  expect_equal(unname(fill$palette.cache), maf_colours(4))

  maf <- c(0, 0.0005, 0.001, 0.0099, 0.01, 0.05, 0.3, 0.5, NA, -0.1)
  cls <- .maf_classes(maf, c(0.001, 0.01, 0.05))
  expect_identical(as.character(cls),
                   c("< 0.1%", "< 0.1%", "0.1-1%", "0.1-1%", "1-5%", ">= 5%", ">= 5%",
                     ">= 5%", NA, NA))
  expect_identical(levels(cls), c(">= 5%", "1-5%", "0.1-1%", "< 0.1%"))
})

test_that("frequencies can be coloured on a continuous scale", {
  p <- manhattan(gw, maf_scale = "continuous")
  built <- ggplot2::ggplot_build(p[[1]])
  fill <- built$plot$scales$get_scales("fill")
  expect_s3_class(fill, "ScaleContinuous")
  expect_equal(10^fill$get_limits()[2], 0.1)
  pts <- point_layer(p[[1]], "solid")
  expect_gt(length(unique(pts$fill)), 3)
  expect_s3_class(manhattan(gw, maf_scale = "continuous", colours = maf_ramp("blue")),
                  "patchwork")
})

test_that("rare variants are drawn on top and unknown frequencies underneath", {
  na <- gw
  na$EAF[na$P < 5e-8][1:5] <- NA
  p <- manhattan(na)
  pts <- point_layer(p[[1]], "solid")
  pal <- c(maf_colours(4), .ink$unknown)
  rank <- match(pts$fill, pal)                      # 1 common ... 4 ultra-rare, 5 unknown
  rank[rank == 5] <- 0
  expect_false(is.unsorted(rank))
})

test_that("a GWAS with no hits still plots", {
  null <- gw
  null$P <- runif(nrow(null), 1e-4, 1)
  expect_s3_class(manhattan(null), "ggplot")
})

test_that("odds ratios and linear effects plot", {
  or <- gw
  or$OR <- exp(or$BETA)
  or$BETA <- NULL
  expect_s3_class(manhattan(or, effect_scale = "linear"), "patchwork")
})

test_that("missing frequencies become an explicit class", {
  na <- gw
  na$EAF[na$P < 5e-8][1:3] <- NA
  p <- manhattan(na)
  fill <- ggplot2::ggplot_build(p[[1]])$plot$scales$get_scales("fill")
  expect_true("Unknown" %in% fill$get_limits())
})

test_that("simulate_gwas is reproducible and leaves the RNG alone", {
  set.seed(99)
  before <- runif(1)
  set.seed(99)
  a <- simulate_gwas(n_variants = 1e4, seed = 3)
  after <- runif(1)
  expect_equal(before, after)
  b <- simulate_gwas(n_variants = 1e4, seed = 3)
  expect_identical(a, b)
  expect_named(a, c("CHR", "BP", "SNP", "A1", "A2", "EAF", "BETA", "SE", "P"))
})
