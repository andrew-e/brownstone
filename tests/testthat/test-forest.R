gw <- simulate_gwas(n_variants = 5e4, seed = 1)

forest_points <- function(p) {
  i <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1)))
  ggplot2::layer_data(p, i)
}

test_that("the forest plot builds and renders on the default pdf() device", {
  p <- forest_plot(gw)
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit({ grDevices::dev.off(); unlink(f) }, add = TRUE)
  expect_no_error(print(p))
})

test_that("one row per lead variant, largest effect at the top", {
  loci <- find_loci(gw)
  p <- forest_plot(gw, effect_scale = "linear")
  pts <- forest_points(p)
  expect_equal(nrow(pts), nrow(loci))
  expect_equal(sort(pts$x), sort(loci$beta))
  top_to_bottom <- pts$x[order(-pts$y)]
  expect_false(is.unsorted(rev(abs(top_to_bottom))))
  # 95% confidence intervals around each estimate.
  i <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomLinerange"), logical(1)))
  ci <- ggplot2::layer_data(p, i)
  ci <- ci[order(ci$y), ]
  pts <- pts[order(pts$y), ]
  expect_equal((ci$xmin + ci$xmax) / 2, pts$x)
  expect_equal((ci$xmax - ci$xmin) / (2 * 1.96), loci$se[match(pts$x, loci$beta)])
})

test_that("rows are labelled with the nearest gene, and its variant ID on the right", {
  row_labels <- function(p, side = "y") {
    vapply(ggplot2::ggplot_build(p)$layout$panel_params[[1]][[side]]$get_labels(),
           function(e) paste(deparse(e), collapse = ""), "")
  }
  loci <- find_loci(gw)
  p <- forest_plot(gw)
  lab <- row_labels(p)
  top <- loci[which.max(abs(loci$beta)), ]
  expect_equal(lab[length(lab)], paste0("italic(\"", top$gene, "\")"))
  ids <- row_labels(p, "y.sec")
  expect_equal(ids[length(ids)], paste0("\"", .short_id(top$snp), "\""))
  # The IDs are muted, the genes in ink.
  expect_equal(p$theme$axis.text.y.right$colour, .ink$muted)
  expect_equal(p$theme$axis.text.y$colour, .ink$primary)
  # Without genes the ID is the row's only label, on the left.
  bare <- forest_plot(gw, annotate = NULL)
  expect_false(any(grepl("italic", row_labels(bare))))
  expect_s3_class(.forest_id_axis(list(right = c("''", "''"))), "waiver")
})

test_that("rows are capped at the most significant loci", {
  expect_message(p <- forest_plot(gw, max_loci = 5, effect_scale = "linear"),
                 "Showing the 5 most significant")
  loci <- find_loci(gw, annotate = FALSE)
  expect_setequal(forest_points(p)$x, loci$beta[1:5])
})

test_that("odds ratios are drawn on the log scale and labelled as odds ratios", {
  or <- gw
  or$OR <- exp(or$BETA)
  or$BETA <- NULL
  p <- forest_plot(or, effect_scale = "linear")
  expect_equal(p$labels$x, "OR (minor allele)")
  expect_equal(sort(forest_points(p)$x), sort(find_loci(or, annotate = FALSE)$beta))
})

test_that("no loci gives an empty plot; no effects, an error", {
  null <- gw
  set.seed(5)
  null$P <- stats::runif(nrow(null), 1e-4, 1)
  expect_message(p <- forest_plot(null), "empty")
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
  expect_error(suppressWarnings(forest_plot(gw[, c("CHR", "BP", "P")])), "effect size")
})

test_that("chr:pos:ref:alt IDs are shortened beside gene names", {
  expect_equal(.short_id(c("11:116828494:A:G", "chr1_55039974_G_T", "rs123", "X:100:AT:A")),
               c("11:116828494", "chr1_55039974", "rs123", "X:100"))
})

# forest_plot(compare = )

gw2 <- replicate_gwas(gw)

compare_points <- function(p, which = c("solid", "hollow")) {
  which <- match.arg(which)
  param <- if (which == "solid") "colour" else "fill"
  i <- which(vapply(p$layers, function(l) {
    inherits(l$geom, "GeomPoint") && identical(l$aes_params[[param]], "white")
  }, logical(1)))
  ggplot2::layer_data(p, i)
}
# Estimates by row: the first GWAS is drawn just above each row, the second
# just below.
estimates <- function(p) {
  pts <- rbind(compare_points(p, "solid")[, c("x", "y")],
               compare_points(p, "hollow")[, c("x", "y")])
  data.frame(row = round(pts$y), gwas = ifelse(pts$y > round(pts$y), 1L, 2L), beta = pts$x)
}
pooled_regions <- function(a, b) {
  cols <- c("chr", "pos", "snp", "logp")
  pool <- rbind(find_loci(a, annotate = FALSE)[, cols], find_loci(b, annotate = FALSE)[, cols])
  pool[.clump(pool, -log10(5e-8), 5e5)$lead, ]
}

test_that("a comparison builds and renders on the default pdf() device", {
  p <- forest_plot(gw, compare = gw2, effect_scale = "linear")
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit({ grDevices::dev.off(); unlink(f) }, add = TRUE)
  expect_no_error(print(p))
  expect_equal(p$scales$get_scales("shape")$get_limits(), c("gw", "gw2"))
})

test_that("one row per region, with both GWAS's estimates of the same variant", {
  reps <- pooled_regions(gw, gw2)
  est <- estimates(forest_plot(gw, compare = gw2, effect_scale = "linear"))
  expect_equal(length(unique(est$row)), nrow(reps))
  expect_true(all(table(est$row) == 2))
  d1 <- standardise_sumstats(gw)
  d2 <- standardise_sumstats(gw2)
  i <- match(reps$snp, d1$snp)
  expect_setequal(round(est$beta[est$gwas == 1], 8), round(d1$beta[i], 8))
  expect_setequal(round(est$beta[est$gwas == 2], 8), round(d2$beta[i], 8))
  # Rows run from the largest |beta| (top) down.
  big <- tapply(abs(est$beta), est$row, max)
  expect_false(is.unsorted(big))
})

test_that("estimates that don't reach the threshold in their GWAS are hollow", {
  p <- forest_plot(gw, compare = gw2, effect_scale = "linear")
  hollow <- compare_points(p, "hollow")
  expect_gt(nrow(hollow), 0)
  d2 <- standardise_sumstats(gw2)
  weak <- d2$beta[d2$logp < -log10(5e-8)]
  expect_true(all(round(hollow$x, 8) %in% round(weak, 8)))
})

test_that("the second GWAS is aligned to the first's alleles", {
  # The same GWAS reported against the other allele.
  swapped <- gw2
  swapped$A1 <- gw2$A2
  swapped$A2 <- gw2$A1
  swapped$BETA <- -gw2$BETA
  swapped$EAF <- 1 - gw2$EAF
  as_is <- function(b) {
    est <- estimates(forest_plot(gw, compare = b, effect_allele = "as_is",
                                 effect_scale = "linear"))
    est[order(est$row, est$gwas), ]
  }
  expect_equal(as_is(swapped), as_is(gw2))

  # A strand flip at one variant: its second estimate is left out.
  reps <- pooled_regions(gw, gw2)
  flip <- c(A = "T", C = "G", G = "C", T = "A")
  k <- gw2$SNP == reps$snp[1]
  strand <- gw2
  strand$A1[k] <- flip[gw2$A1[k]]
  strand$A2[k] <- flip[gw2$A2[k]]
  msg <- capture_messages(p <- forest_plot(gw, compare = strand, effect_scale = "linear"))
  expect_match(msg, "Alleles of 1 variant", all = FALSE)
  expect_equal(sum(table(estimates(p)$row) == 1), 1)
})

test_that("without alleles, common variants whose minor allele may differ are listed", {
  expect_message(forest_plot(gw[, names(gw) != "A2"], compare = gw2[, names(gw2) != "A2"]),
                 "may differ")
})

test_that("a variant missing from one GWAS gets a single point", {
  reps <- pooled_regions(gw, gw2)
  expect_message(p <- forest_plot(gw, compare = gw2[gw2$SNP != reps$snp[2], ],
                                  effect_scale = "linear"),
                 "1 of [0-9]+ variants have no estimate in")
  expect_equal(sum(table(estimates(p)$row) == 1), 1)
})

test_that("odds ratios can't be compared with betas", {
  or <- gw2
  or$OR <- exp(or$BETA)
  or$BETA <- NULL
  expect_error(forest_plot(gw, compare = or), "odds ratios")
  expect_error(forest_plot(gw, compare = gw2, names = "one"), "names")
})
