gw <- simulate_gwas(n_variants = 5e4, seed = 1)

point_data <- function(p) {
  i <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1)))
  ggplot2::layer_data(p, i)
}

legend_text <- function(p) {
  labels <- ggplot2::ggplot_build(p)$plot$scales$get_scales("colour")$get_labels()
  vapply(labels, function(e) paste(deparse(e), collapse = ""), "")
}

test_that("the QQ plot builds and renders on the default pdf() device", {
  p <- qq_plot(gw)
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit({ grDevices::dev.off(); unlink(f) }, add = TRUE)
  expect_no_error(print(p))
})

test_that("genomic inflation is the median chi-square over its expectation", {
  expect_equal(.lambda_gc(-log10(0.5)), 1)
  # A null GWAS simulated with inflation 1.10, overall and in every class.
  null <- simulate_gwas(n_variants = 2e5, loci = sim_loci[0, ], n_background = 0,
                        inflation = 1.10, seed = 1)
  d <- standardise_sumstats(null)
  expect_equal(.lambda_gc(stats::median(d$logp)), 1.10, tolerance = 0.02 / 1.10)
  cls <- .maf_classes(d$maf, c(0.001, 0.01, 0.05))
  by_class <- tapply(d$logp, cls, function(x) .lambda_gc(stats::median(x)))
  expect_true(all(abs(by_class - 1.10) < 0.03))
})

test_that("lambda is printed on the plot, overall and per class", {
  null <- simulate_gwas(n_variants = 2e5, loci = sim_loci[0, ], n_background = 0,
                        inflation = 1.10, seed = 1)
  p <- qq_plot(null)
  text <- Filter(function(l) inherits(l$geom, "GeomText"), p$layers)
  expect_length(text, 1)
  overall <- .lambda_gc(stats::median(standardise_sumstats(null)$logp))
  expect_equal(text[[1]]$aes_params$label,
               sprintf("lambda[GC] == '%s'", formatC(overall, format = "f", digits = 3)))
  expect_true(all(grepl("lambda", legend_text(p))))
  # Classes too small for a precise lambda get no number.
  tiny <- null
  tiny$EAF[seq_len(nrow(tiny)) > 500 & pmin(tiny$EAF, 1 - tiny$EAF) < 0.001] <- 0.3
  expect_false(grepl("lambda", legend_text(qq_plot(tiny))[4]))
  quiet <- qq_plot(null, lambda = FALSE)
  expect_length(Filter(function(l) inherits(l$geom, "GeomText"), quiet$layers), 0)
})

test_that("each frequency class is compared with its own expectation", {
  p <- qq_plot(gw)
  colour <- ggplot2::ggplot_build(p)$plot$scales$get_scales("colour")
  expect_equal(colour$get_limits(), c(">= 5%", "1-5%", "0.1-1%", "< 0.1%"))
  pts <- point_data(p)
  # The most significant variant of each class is at -log10(0.5 / n_class).
  d <- standardise_sumstats(gw)
  n_class <- table(.maf_classes(d$maf, c(0.001, 0.01, 0.05)))
  expect_equal(sort(as.numeric(tapply(pts$x, pts$colour, max))),
               sort(-log10(0.5 / as.numeric(n_class))))
  # Thinned, but every point beyond the significance line is drawn.
  expect_lt(nrow(pts), nrow(gw))
  expect_equal(sum(pts$y >= -log10(5e-8) & pts$y < 20), sum(d$logp >= -log10(5e-8) & d$logp < 20))
  # The y axis breaks at 20, as in manhattan().
  expect_true(any(vapply(p$layers, function(l) identical(l$aes_params$fill, "white"), logical(1))))

  one <- qq_plot(gw, by = "none")
  expect_equal(max(point_data(one)$x), -log10(0.5 / nrow(gw)))
  expect_length(unique(point_data(one)$colour), 1)
})

test_that("missing frequencies form their own class; no frequencies, one curve", {
  na <- gw
  na$EAF[1:3000] <- NA
  colour <- ggplot2::ggplot_build(qq_plot(na))$plot$scales$get_scales("colour")
  expect_true("Unknown" %in% colour$get_limits())
  expect_message(p <- qq_plot(gw[, c("CHR", "BP", "P")]), "one curve")
  expect_length(unique(point_data(p)$colour), 1)
})

test_that("filtered or degenerate P-values are flagged", {
  filtered <- gw[gw$P < 0.2, ]
  expect_warning(p <- qq_plot(filtered), "look filtered")
  expect_length(Filter(function(l) inherits(l$geom, "GeomText"), p$layers), 0)
  ones <- gw
  ones$P[1:1000] <- 1
  expect_warning(qq_plot(ones), "exactly 1")
  expect_silent(qq_plot(gw))
})

test_that("simulate_gwas() can simulate a GWAS with no causal loci", {
  null <- simulate_gwas(n_variants = 1e4, loci = sim_loci[0, ], n_background = 0, seed = 2)
  expect_equal(nrow(attr(null, "causal")), 0)
  expect_lt(min(null$P[-which.min(null$P)]), 1)
})

test_that("any column can define the strata", {
  set.seed(3)
  gw$INFO_BIN <- factor(sample(c("< 0.3", "0.3-0.8", ">= 0.8"), nrow(gw), TRUE),
                        levels = c("< 0.3", "0.3-0.8", ">= 0.8"))
  gw$INFO_BIN[1:100] <- NA
  p <- qq_plot(gw, by = "INFO_BIN", keep = "INFO_BIN")
  colour <- ggplot2::ggplot_build(p)$plot$scales$get_scales("colour")
  expect_equal(colour$get_limits(), c("< 0.3", "0.3-0.8", ">= 0.8", "Unknown"))
  expect_equal(colour$name, "INFO_BIN")
  expect_equal(unname(colour$palette.cache[1:3]), maf_colours(3))
  # Each group is compared with its own expectation.
  pts <- point_data(p)
  n_group <- as.numeric(table(gw$INFO_BIN, useNA = "always"))
  expect_equal(sort(as.numeric(tapply(pts$x, pts$colour, max))),
               sort(-log10(0.5 / n_group)))
})

test_that("more than four strata get a panel each, with their lambda", {
  p <- qq_plot(gw, by = "chr")
  expect_s3_class(p$facet, "FacetWrap")
  strips <- levels(ggplot2::layer_data(p, 1)$PANEL)
  expect_length(strips, 22)
  strip <- as.character(ggplot2::ggplot_build(p)$layout$layout$strip)
  # Panels follow the karyotype, not alphabetical order, each with its lambda.
  expect_match(strip[1], "^'1' ~~~ lambda == '1\\.[0-9]{3}'$")
  expect_match(strip[10], "^'10' ~~~ lambda")
  # ...where a chromosome has the 1,000 variants needed for a precise one.
  expect_equal(grepl("lambda", strip), as.vector(table(gw$CHR) >= 1000))
  expect_match(paste(deparse(p$labels$subtitle), collapse = ""), "lambda[GC]", fixed = TRUE)
})

test_that("an unknown or continuous `by` column is refused", {
  expect_error(qq_plot(gw, by = "INFO"), "keep = \"INFO\"")
  expect_error(qq_plot(gw, by = "pos"), "bin a continuous column")
})
