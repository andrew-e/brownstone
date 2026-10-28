gw <- simulate_gwas(n_variants = 2e5, seed = 1)
sim_ld <- attr(gw, "ld")
causal <- attr(gw, "causal")
at <- function(g) {
  x <- causal[causal$gene == g, ]
  paste0(x$chr, ":", x$pos)
}
x_range <- function(p) ggplot2::ggplot_build(p[[1]])$layout$panel_params[[1]]$x.range
# The point layers of a panel (points, or ggstar shapes for frequency
# classes), by their fixed parameters.
is_point <- function(l) inherits(l$geom, c("GeomPoint", "GeomStar"))
lead_layer <- function(panel) {
  i <- which(vapply(panel$layers, function(l) identical(l$aes_params$colour, .ink$primary) &&
                      is_point(l), logical(1)))
  ggplot2::layer_data(panel, i)
}
point_rows <- function(panel) {
  sum(vapply(seq_along(panel$layers), function(i) {
    if (is_point(panel$layers[[i]])) nrow(ggplot2::layer_data(panel, i)) else 0L
  }, 0L))
}

test_that("a regional plot builds and renders, with LD, effect and gene tracks", {
  expect_message(p <- region_plot(gw, snp = at("HMGCR"), ld = sim_ld), "your table")
  expect_s3_class(p, "patchwork")
  for (k in 1:3) expect_silent(ggplot2::ggplot_build(p[[k]]))
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f, width = 8, height = 7)
  on.exit({ grDevices::dev.off(); unlink(f) }, add = TRUE)
  expect_no_error(print(p))
})

test_that("the region is a variant, a gene or a range, or the top hit", {
  d <- standardise_sumstats(gw)
  h <- causal[causal$gene == "HMGCR", ]
  quiet <- function(...) suppressMessages(region_plot(gw, ..., use_ld_info = FALSE))
  expect_equal(x_range(quiet(snp = at("HMGCR"))), h$pos + c(-2.5e5, 2.5e5))
  expect_equal(x_range(quiet(snp = at("HMGCR"), flank = 1e5)), h$pos + c(-1e5, 1e5))
  g <- genes_grch38[genes_grch38$gene == "HMGCR", ]
  expect_equal(x_range(quiet(gene = "HMGCR")), c(g$start - 2.5e5, g$end + 2.5e5))
  expect_equal(x_range(quiet(chr = "chr5", start = 75e6, end = 75.5e6)), c(75e6, 75.5e6))
  top <- d$pos[which.max(d$logp)]
  expect_equal(x_range(quiet()), top + c(-2.5e5, 2.5e5))

  expect_error(region_plot(gw, snp = at("HMGCR"), gene = "HMGCR"), "one of")
  expect_error(region_plot(gw, gene = "NOTAGENE", use_ld_info = FALSE), "not in the GRCh38 gene table")
  expect_error(region_plot(gw, snp = "rs999", use_ld_info = FALSE), "not in the data")
  expect_error(region_plot(gw, chr = 5, start = 1e6, use_ld_info = FALSE), "together")
})

test_that("r2 is binned, with frequency-capped variants in the lowest bin", {
  cls <- .ld_classes(c(1, 0.85, 0.5, 0.1, NA, NA), maf = c(0.3, 0.3, 0.3, 0.3, 0.001, 0.3),
                     lead_maf = 0.3, breaks = c(0.2, 0.4, 0.6, 0.8), has_maf = TRUE)
  expect_equal(as.character(cls), c(">= 0.8", ">= 0.8", "0.4-0.6", "< 0.2", "< 0.2", "No LD data"))
  expect_equal(levels(cls), c("< 0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", ">= 0.8", "No LD data"))

  p <- suppressMessages(region_plot(gw, snp = at("HMGCR"), ld = sim_ld))
  # The lead is in the darkest bin, in both panels.
  expect_equal(lead_layer(p[[1]])$fill, .ld_colours()[5])
  expect_equal(lead_layer(p[[2]])$fill, .ld_colours()[5])
  fill <- ggplot2::ggplot_build(p[[1]])$plot$scales$get_scales("fill")
  expect_equal(fill$get_limits(), levels(cls))
})

test_that("LD from LDlink is matched to the region, and a build mismatch is caught", {
  d <- standardise_sumstats(gw)
  h <- causal[causal$gene == "HMGCR", ]
  rows <- sim_ld[sim_ld$BP_A == h$pos, ]
  tab <- data.frame(chr = as.character(rows$CHR_B), pos = rows$BP_B, snp = NA_character_,
                    alleles = NA, maf = NA, distance = rows$BP_B - h$pos, dprime = NA,
                    r2 = rows$R2)
  local_mocked_bindings(ld_proxy = function(snp, ...) tab)
  expect_message(region_plot(gw, snp = at("HMGCR"), token = "x"),
                 paste0("r² for ", sum(abs(tab$distance) <= 2.5e5) + 1L))

  # Queried by rsID, LDlink places the lead elsewhere: positions can't be trusted.
  rs <- gw
  rs$SNP <- paste0("rs", seq_len(nrow(rs)))
  lead <- rs$SNP[rs$CHR == as.integer(h$chr) & rs$BP == h$pos]
  moved <- rbind(data.frame(chr = h$chr, pos = h$pos + 1000, snp = lead, alleles = NA, maf = NA,
                            distance = 0, dprime = 1, r2 = 1), tab)
  local_mocked_bindings(ld_proxy = function(snp, ...) moved)
  expect_warning(suppressMessages(region_plot(rs, snp = lead, token = "x")), "is `build` right")
})

test_that("without LD the plot is a zoomed-in Manhattan plot", {
  expect_error(region_plot(gw, snp = at("HMGCR"), token = ""), "use_ld_info = FALSE")
  expect_error(region_plot(gw, snp = at("HMGCR"), ld = NULL), "use_ld_info = FALSE")
  p <- suppressMessages(region_plot(gw, snp = at("HMGCR"), token = "", use_ld_info = FALSE))
  solid <- Filter(function(l) inherits(l$geom, "GeomPoint") &&
                    identical(l$aes_params$colour, "white"), p[[1]]$layers)
  expect_gt(length(solid), 0)
  fill <- ggplot2::ggplot_build(p[[1]])$plot$scales$get_scales("fill")
  expect_equal(fill$name, "Minor allele frequency")
  expect_message(region_plot(gw, snp = at("HMGCR"), ld = data.frame(SNP_A = "x", SNP_B = "y",
                                                                     R2 = 1)), "no rows for")
})

test_that("the effect track shows every variant above the suggestive line", {
  d <- standardise_sumstats(gw)
  h <- causal[causal$gene == "HMGCR", ]
  reg <- d[d$chr == h$chr & abs(d$pos - h$pos) <= 2.5e5, ]
  p <- suppressMessages(region_plot(gw, snp = at("HMGCR"), ld = sim_ld))
  n_hits <- sum(reg$logp >= 5 & is.finite(reg$beta) & is.finite(reg$se))
  expect_equal(point_rows(p[[2]]), n_hits)   # the lead is one of them here
  p2 <- suppressMessages(region_plot(gw, snp = at("HMGCR"), ld = sim_ld, suggestive = 1e-3))
  expect_equal(point_rows(p2[[2]]), sum(reg$logp >= 3 & is.finite(reg$beta)))
})

test_that("an LD column of the data can be used, and a rare lead is flagged", {
  gw$R2 <- runif(nrow(gw))
  expect_s3_class(suppressMessages(region_plot(gw, snp = at("HMGCR"), ld = "R2")), "patchwork")
  rare <- causal[causal$gene == "APOC3", ]
  near <- gw$BP[gw$CHR == as.integer(rare$chr) & abs(gw$BP - rare$pos) < 2e5][1:2]
  tab <- data.frame(chr = rare$chr, pos = near, r2 = 0.5)
  expect_message(region_plot(gw, snp = at("APOC3"), ld = tab, genes = FALSE), "is rare")
})

test_that("frequency classes are circle, square, triangle, star, rarest most distinct", {
  lev <- c(">= 5%", "1-5%", "0.1-1%", "< 0.1%", "Unknown")
  expect_equal(unname(.maf_shapes(lev)), c(15, 13, 11, 1, 29))   # ggstar starshapes
  expect_error(.maf_shapes(c("a", "b", "c", "d", "e")), "at most four")
  p <- suppressMessages(region_plot(gw, snp = at("HMGCR"), ld = sim_ld))
  shape <- ggplot2::ggplot_build(p[[1]])$plot$scales$get_scales("starshape")
  expect_equal(shape$name, "Minor allele frequency")
  # With only common hits every shape is the plain circle: no key for it.
  lev <- c(">= 5%", "1-5%", "0.1-1%", "< 0.1%")
  r <- data.frame(class = factor(lev[c(1, 1, 3)], levels = lev), hit = c(TRUE, TRUE, FALSE),
                  lead = c(TRUE, FALSE, FALSE))
  ldc <- factor(c("a", "b", "c"), levels = c("a", "b", "No LD data"))
  guide <- function(r) .region_scales(ldc, r, legend = TRUE)[[2]]$guide
  expect_identical(guide(r), "none")
  r$hit[3] <- TRUE
  expect_s3_class(guide(r), "Guide")
})

test_that("the -log10(P) panel stands on a ground line, clear of the track", {
  p <- suppressMessages(region_plot(gw, snp = at("HMGCR"), ld = sim_ld))
  expect_s3_class(p[[1]]$theme$axis.line.x.bottom, "element_line")
  y <- ggplot2::ggplot_build(p[[1]])$layout$panel_params[[1]]$y.range
  expect_lt(y[1], 0)   # the floor stands above the line
  b <- manhattan(gw)
  expect_s3_class(b[[1]]$theme$axis.line.x.bottom, "element_line")
  expect_lt(ggplot2::ggplot_build(b[[1]])$layout$panel_params[[1]]$y.range[1], 0)
  hanging <- manhattan(gw, track_direction = "down")
  expect_equal(ggplot2::ggplot_build(hanging[[1]])$layout$panel_params[[1]]$y.range[1], 0)
})

test_that("genes are packed into rows, nearest the lead first", {
  # A and B overlap; C is far enough away for its name not to collide with A's.
  g <- data.frame(gene = c("A", "B", "C"), start = c(100, 150, 5000), end = c(300, 350, 6000),
                  strand = "+")
  pk <- .pack_genes(g, start = 0, end = 1e5, lead_pos = 200)
  expect_equal(pk$row, c(1L, 2L, 1L))
  expect_equal(.pack_genes(g, start = 0, end = 1e5, lead_pos = 200, max_rows = 1)$row,
               c(1L, NA, 1L))
  expect_equal(nrow(.pack_genes(g[0, ], 0, 1e5, 200)), 0)
  # A name that would cross the window's edge is aligned to it.
  edge <- .pack_genes(g[1, ], start = 0, end = 1e5, lead_pos = 200)
  expect_equal(c(edge$label_x, edge$hjust), c(0, 0))

  p <- suppressMessages(region_plot(gw, chr = 5, start = 140.7e6, end = 141.3e6,
                                    use_ld_info = FALSE))
  txt <- unlist(lapply(p[[3]]$layers, function(l) l$aes_params$label))
  expect_true(any(grepl("^\\+[0-9]+ genes$", txt)))
  g13 <- genes_grch38[genes_grch38$chr == "13", ]
  edges <- sort(c(g13$start, g13$end))
  desert <- edges[which.max(diff(edges))]
  p <- suppressMessages(region_plot(gw, chr = 13, start = desert + 1e5, end = desert + 6e5,
                                    use_ld_info = FALSE))
  txt <- unlist(lapply(p[[3]]$layers, function(l) l$aes_params$label))
  expect_true("No protein-coding genes" %in% txt)
  expect_error(suppressMessages(region_plot(gw, snp = at("HMGCR"), use_ld_info = FALSE,
                                            genes = FALSE))[[3]])
})

test_that("simulated GWAS carry the LD they were built with", {
  expect_named(sim_ld, c("CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2"))
  expect_true(all(sim_ld$R2 > 0.01 & sim_ld$R2 <= 1))
  expect_true(all(sim_ld$SNP_A %in% gw$SNP) && all(sim_ld$SNP_B %in% gw$SNP))
})
