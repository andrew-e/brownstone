#' Forest plot of the lead variant of each locus
#'
#' @description
#' One row per genome-wide significant locus (found by [find_loci()]): the lead
#' variant's effect with its 95% confidence interval, largest effects at the
#' top.
#'
#' * **Effects are signed**, oriented to the minor allele where an allele
#'   frequency is available (see `effect_allele` in [standardise_sumstats()]).
#'   Odds ratios are drawn on the log scale, so that intervals are symmetric,
#'   and labelled as odds ratios.
#' * **Lead variants are coloured by minor allele frequency**, as in
#'   [manhattan()].
#' * **Rows are labelled** with the nearest gene and the variant ID, or the ID
#'   alone with `annotate = NULL`.
#' * **The effect axis reacts to the data**: linear, unless a few effects dwarf
#'   the rest; then, as on the effect track of [manhattan()], it breaks and puts
#'   the outliers on a log scale. Confidence intervals beyond the axis are
#'   clipped at its edge.
#'
#' Lead-variant effects are not corrected for winner's curse, which inflates
#' rare variants and loci just past the threshold the most.
#'
#' @section Comparing two GWAS:
#' With `compare`, each row is one variant with both GWAS's estimates side by
#' side: the first GWAS above (circles), the second below (diamonds), and
#' alternate variants shaded so that each pair reads as a group.
#'
#' * **Which variants.** Loci are found in each GWAS, and loci of the two within
#'   `window` of each other are merged: each region is represented by the lead
#'   variant of whichever GWAS is more significant there, and that same variant
#'   is looked up in the other GWAS (by ID, or by position when it is the only
#'   variant there).
#' * **Hollow points** are estimates that don't reach `threshold` in their own
#'   GWAS, so a signal found in one GWAS and only nominal in the other stands
#'   out. A variant missing from one GWAS has a single point.
#' * **Alleles.** Effects are compared per the same allele. When both GWAS have
#'   effect and other allele columns (detected, or named with `ea =` and
#'   `oa =`), the second GWAS's effects are aligned to the first's alleles, and
#'   variants whose alleles don't match (a strand flip, or a different variant)
#'   lose the second estimate. Without alleles, both are per the minor allele
#'   of each GWAS (with allele frequencies) or as reported; variants with a
#'   minor allele frequency above 40% in either are listed, since their minor
#'   alleles may differ.
#'
#' @inheritParams find_loci
#' @inheritParams volcano_plot
#' @param max_loci Show at most this many loci (the most significant), still
#'   ordered by effect size. Each row needs about 0.15 inches of figure height
#'   (about 0.25 with `compare`).
#' @param effect_scale How the effect axis reacts to the data: `"auto"` (the
#'   default; linear, breaking for outlying effects) or `"linear"`.
#' @param compare Optional summary statistics of a second GWAS, drawn beside
#'   the first for every locus of either (see "Comparing two GWAS"). Column
#'   names in `...` apply to both; when the files differ, standardise each with
#'   [standardise_sumstats()] first (naming `ea` and `oa` to keep alleles).
#' @param names Names of the two GWAS for the legend, with `compare`. By
#'   default, the expressions passed as `data` and `compare`.
#'
#' @return A ggplot.
#' @export
#' @examples
#' gwas <- simulate_gwas(n_variants = 2e5)
#' forest_plot(gwas)
#'
#' # A second GWAS of the same variants: weaker signals, its own noise.
#' replication <- gwas
#' z <- 0.75 * gwas$BETA / gwas$SE + sqrt(1 - 0.75^2) * rnorm(nrow(gwas))
#' replication$BETA <- z * gwas$SE
#' replication$P <- 2 * pnorm(-abs(z))
#' forest_plot(gwas, compare = replication, max_loci = 20)
forest_plot <- function(data, ..., compare = NULL, names = NULL,
                        build = c("GRCh38", "GRCh37"), threshold = 5e-8,
                        window = 5e5, annotate = "nearest_gene", max_loci = 40,
                        effect_scale = c("auto", "linear"),
                        maf_breaks = c(0.001, 0.01, 0.05), colours = NULL,
                        base_size = 10, point_size = 1) {
  build <- match.arg(build)
  effect_scale <- match.arg(effect_scale)
  if (!is.null(compare)) {
    names <- names %||% c(deparse1(substitute(data)), deparse1(substitute(compare)))
    if (length(names) != 2L) stop("`names` must name the two GWAS.", call. = FALSE)
    return(.forest_compare(data, compare, ..., names = as.character(names), build = build,
                           threshold = threshold, window = window, annotate = annotate,
                           max_loci = max_loci, effect_scale = effect_scale,
                           maf_breaks = maf_breaks, colours = colours,
                           base_size = base_size, point_size = point_size))
  }
  d <- .standardise(data, ..., ids = FALSE)
  if (!isTRUE(attr(d, "has_effect"))) {
    stop("forest_plot() needs effect size and standard error columns; ",
         "name them with `beta =` (or `or =`) and `se =`.", call. = FALSE)
  }
  log_or <- identical(attr(d, "effect_scale"), "log_or")
  allele <- attr(d, "effect_allele")
  has_maf <- isTRUE(attr(d, "has_maf"))
  sig <- -log10(threshold)

  loci <- find_loci(d, threshold = threshold, window = window, annotate = annotate,
                    build = build)
  loci <- loci[is.finite(loci$beta) & is.finite(loci$se), ]
  if (!nrow(loci)) return(.forest_empty(threshold, log_or, allele, base_size))
  if (nrow(loci) > max_loci) {
    message("Showing the ", max_loci, " most significant of ", nrow(loci),
            " loci; set `max_loci` to show more.")
    loci <- loci[seq_len(max_loci), ]   # find_loci() lists them by significance
  }

  # Rows from the largest |beta| (top) down; factor levels run bottom to top.
  loci <- loci[order(abs(loci$beta)), ]
  gene <- if ("gene" %in% names(loci)) loci$gene else rep(NA_character_, nrow(loci))
  row_lab <- .forest_row_labels(gene, loci$snp)
  loci$row <- factor(seq_len(nrow(loci)))
  loci$lo <- loci$beta - 1.96 * loci$se
  loci$hi <- loci$beta + 1.96 * loci$se

  cp <- .class_palette(loci, has_maf, maf_breaks, colours, floor_logp = sig)
  loci$class <- cp$class

  # Room for every estimate, but one very wide interval doesn't squash the rest
  # (as on manhattan()'s effect track).
  top <- max(min(max(abs(c(loci$lo, loci$hi))), 1.6 * max(abs(loci$beta))), 1e-12)
  ax <- .effect_axis(abs(loci$beta), top, effect_scale, max_estimate = max(abs(loci$beta)))
  pts <- point_size * base_size / 10

  ggplot2::ggplot(loci, ggplot2::aes(y = .data$row)) +
    .effect_break(ax, top, signed = TRUE, axis = "x") +
    ggplot2::geom_vline(xintercept = 0, colour = .ink$axis, linewidth = 0.4) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = .data$lo, xmax = .data$hi),
                            colour = .ink$muted, linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(x = .data$beta, fill = .data$class), shape = 21,
                        colour = "white", stroke = 0.3 * pts, size = 2.4 * pts) +
    ggplot2::scale_fill_manual(
      values = cp$pal, labels = .maf_legend_labels(names(cp$pal)), drop = FALSE,
      name = "Minor allele frequency",
      guide = if (has_maf) {
        ggplot2::guide_legend(override.aes = list(size = 3))
      } else "none"
    ) +
    ggplot2::scale_x_continuous(
      transform = ax$trans, limits = c(-top, top), oob = scales::oob_squish,
      breaks = .track_breaks(ax, top, signed = TRUE, log_or = log_or),
      labels = function(b) .as_plotmath(.effect_labels(b, log_or)),
      expand = ggplot2::expansion(mult = 0.04)
    ) +
    ggplot2::scale_y_discrete(labels = .as_plotmath(row_lab$left),
                              sec.axis = .forest_id_axis(row_lab)) +
    ggplot2::labs(x = .effect_title(signed = TRUE, log_or, allele), y = NULL) +
    theme_maffect(base_size = base_size) + .legend_on_top(base_size) +
    .forest_label_theme(base_size) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(colour = .ink$grid, linewidth = 0.3)
    )
}

.forest_empty <- function(threshold, log_or, allele, base_size) {
  p_txt <- parse(text = paste("italic(P) <", .p_plotmath(threshold)))[[1]]
  message("No loci reach P < ", format(threshold), "; the forest plot is empty.")
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, colour = .ink$muted, family = .family(),
                      size = base_size * .text_size$legend / ggplot2::.pt,
                      label = deparse(bquote("No loci reach" ~ .(p_txt))), parse = TRUE) +
    ggplot2::labs(x = .effect_title(signed = TRUE, log_or, allele), y = NULL) +
    theme_maffect(base_size = base_size) +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank())
}

# Row labels, as plotmath: the nearest gene on the left, in ink, and its
# variant ID on the right, in muted ink, so the names read first; a variant
# without a gene has its ID on the left. Beside a gene name, chr:pos:ref:alt
# IDs are shortened to chr:pos: the alleles would make them half the width of
# a single-column figure.
.forest_row_labels <- function(gene, snp) {
  list(left = ifelse(is.na(gene), .pm_string(snp), sprintf("italic(%s)", .pm_string(gene))),
       right = ifelse(is.na(gene), "''", .pm_string(.short_id(snp))))
}

# The right-hand axis of variant IDs, when any row has a gene on the left.
.forest_id_axis <- function(row_lab) {
  if (all(row_lab$right == "''")) return(ggplot2::waiver())
  ggplot2::dup_axis(name = NULL, labels = .as_plotmath(row_lab$right))
}

# Row-label styling shared by the forest plots: genes in ink, IDs muted.
.forest_label_theme <- function(base_size) {
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(colour = .ink$primary),
    axis.text.y.right = ggplot2::element_text(colour = .ink$muted, hjust = 0,
                                              size = base_size * .text_size$note)
  )
}

# forest_plot(compare = ): both GWAS's estimates for one variant per region,
# side by side.
.forest_compare <- function(data, compare, ..., names, build, threshold, window, annotate,
                            max_loci, effect_scale, maf_breaks, colours, base_size,
                            point_size) {
  d <- list(.standardise(data, ..., ids = FALSE, alleles = TRUE),
            .standardise(compare, ..., ids = FALSE, alleles = TRUE))
  for (k in 1:2) {
    if (!isTRUE(attr(d[[k]], "has_effect"))) {
      stop("forest_plot() needs effect size and standard error columns in both GWAS (",
           names[k], " has none); name them with `beta =` (or `or =`) and `se =`.",
           call. = FALSE)
    }
  }
  scales <- vapply(d, function(x) attr(x, "effect_scale"), "")
  if (scales[1] != scales[2]) {
    stop("One GWAS reports odds ratios and the other betas; they can't share an axis.",
         call. = FALSE)
  }
  log_or <- scales[1] == "log_or"
  # With alleles in both, the second GWAS is aligned to the first: effects are
  # per the first GWAS's (minor, if oriented) allele.
  oriented <- vapply(d, function(x) attr(x, "effect_allele"), "")
  aligned <- all(vapply(d, function(x) "allele" %in% names(x), logical(1)))
  allele <- if (aligned) oriented[1] else if (all(oriented == "minor")) "minor" else "as_is"
  has_maf <- any(vapply(d, function(x) isTRUE(attr(x, "has_maf")), logical(1)))
  sig <- -log10(threshold)

  # One variant per region: loci of the two GWAS are pooled and clumped, so
  # each region is represented by the lead of whichever GWAS is stronger there.
  loci <- lapply(d, function(x) {
    l <- find_loci(x, threshold = threshold, window = window, annotate = FALSE)
    l[is.finite(l$beta) & is.finite(l$se), c("chr", "pos", "snp", "logp")]
  })
  pool <- rbind(loci[[1]], loci[[2]])
  if (!nrow(pool)) return(.forest_empty(threshold, log_or, allele, base_size))
  cl <- .clump(pool, sig, window)
  reps <- pool[cl$lead, ]   # most significant first
  if (nrow(reps) > max_loci) {
    message("Showing the ", max_loci, " most significant of ", nrow(reps),
            " loci; set `max_loci` to show more.")
    reps <- reps[seq_len(max_loci), ]
  }
  rownames(reps) <- NULL
  rows <- lapply(d, function(x) .match_variants(reps, x, unique_pos = TRUE))

  # Both GWAS's estimates, long: one row per variant and GWAS.
  est <- do.call(rbind, lapply(1:2, function(k) {
    i <- rows[[k]]
    data.frame(group = seq_len(nrow(reps)), gwas = names[k], k = k,
               beta = d[[k]]$beta[i], se = d[[k]]$se[i], logp = d[[k]]$logp[i],
               maf = d[[k]]$maf[i],
               allele = if ("allele" %in% names(d[[k]])) d[[k]]$allele[i] else NA_character_,
               other = if ("other_allele" %in% names(d[[k]])) d[[k]]$other_allele[i] else NA_character_,
               stringsAsFactors = FALSE)
  }))
  est <- .align_alleles(est, reps, names)
  est <- est[is.finite(est$beta) & is.finite(est$se), ]
  for (k in 1:2) {
    miss <- sum(!seq_len(nrow(reps)) %in% est$group[est$k == k])
    if (miss) message(miss, " of ", nrow(reps), " variants have no estimate in ", names[k], ".")
  }

  # Groups from the largest |beta| (top) down; the first GWAS above the second.
  big <- vapply(seq_len(nrow(reps)), function(g) max(abs(est$beta[est$group == g]), 0), 0)
  rank <- order(order(big))   # 1 = smallest, at the bottom
  est$y <- rank[est$group] + ifelse(est$k == 1, 0.18, -0.18)
  gene <- if (isTRUE(annotate == "nearest_gene")) {
    nearest_gene(reps$chr, reps$pos, build = build)$gene
  } else rep(NA_character_, nrow(reps))
  row_lab <- lapply(.forest_row_labels(gene, reps$snp), function(l) l[order(rank)])
  est$lo <- est$beta - 1.96 * est$se
  est$hi <- est$beta + 1.96 * est$se

  cp <- .class_palette(est, has_maf, maf_breaks, colours, floor_logp = 0)
  est$class <- cp$class
  est$gwas <- factor(est$gwas, levels = names)
  strong <- est[est$logp >= sig, ]
  weak <- est[est$logp < sig, ]

  top <- max(min(max(abs(c(est$lo, est$hi))), 1.6 * max(abs(est$beta))), 1e-12)
  ax <- .effect_axis(abs(strong$beta), top, effect_scale, max_estimate = max(abs(est$beta)))
  pts <- point_size * base_size / 10
  n <- nrow(reps)
  bands <- seq_len(n)[seq_len(n) %% 2 == 0]
  shapes <- stats::setNames(c(21, 23), names)
  p_txt <- parse(text = paste("italic(P) >=", .p_plotmath(threshold)))[[1]]

  ggplot2::ggplot(est, ggplot2::aes(y = .data$y)) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = bands - 0.5, ymax = bands + 0.5,
                      fill = .ink$band, colour = NA) +
    .effect_break(ax, top, signed = TRUE, axis = "x") +
    ggplot2::geom_vline(xintercept = 0, colour = .ink$axis, linewidth = 0.4) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = .data$lo, xmax = .data$hi),
                            colour = .ink$muted, linewidth = 0.6) +
    ggplot2::geom_point(data = weak, ggplot2::aes(x = .data$beta, colour = .data$class,
                                                  shape = .data$gwas),
                        fill = "white", stroke = 0.7 * pts, size = 2 * pts) +
    ggplot2::geom_point(data = strong, ggplot2::aes(x = .data$beta, fill = .data$class,
                                                    shape = .data$gwas),
                        colour = "white", stroke = 0.3 * pts, size = 2.4 * pts,
                        show.legend = TRUE) +
    ggplot2::scale_fill_manual(
      values = cp$pal, labels = .maf_legend_labels(names(cp$pal)), drop = FALSE,
      name = "Minor allele frequency",
      guide = if (has_maf) {
        ggplot2::guide_legend(order = 1, override.aes = list(shape = 21, size = 3))
      } else "none"
    ) +
    ggplot2::scale_colour_manual(values = cp$pal, drop = FALSE, guide = "none") +
    ggplot2::scale_shape_manual(
      values = shapes, limits = names, name = NULL,
      guide = ggplot2::guide_legend(order = 2, override.aes = list(
        fill = .ink$secondary, colour = "white", size = 3))
    ) +
    ggplot2::scale_x_continuous(
      transform = ax$trans, limits = c(-top, top), oob = scales::oob_squish,
      breaks = .track_breaks(ax, top, signed = TRUE, log_or = log_or),
      labels = function(b) .as_plotmath(.effect_labels(b, log_or)),
      expand = ggplot2::expansion(mult = 0.04)
    ) +
    ggplot2::scale_y_continuous(breaks = seq_len(n), labels = .as_plotmath(row_lab$left),
                                limits = c(0.5, n + 0.5), expand = ggplot2::expansion(0),
                                sec.axis = .forest_id_axis(row_lab)) +
    ggplot2::labs(x = .effect_title(signed = TRUE, log_or, allele), y = NULL,
                  caption = bquote("Hollow:" ~ .(p_txt) ~ "in that GWAS")) +
    theme_maffect(base_size = base_size) + .legend_on_top(base_size) +
    .forest_label_theme(base_size) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(colour = .ink$grid, linewidth = 0.3),
      panel.grid.major.y = ggplot2::element_blank(),
      # The two legends' keys on one line, though only the first has a title.
      legend.box.just = "bottom"
    )
}

# Orient the second GWAS's estimates (`est$k == 2`) to the first's alleles,
# group by group. Where either lacks alleles, effects stay as standardised,
# and variants common enough for the minor allele to differ are listed.
.align_alleles <- function(est, reps, names) {
  a <- est[est$k == 1, ]
  b <- est[est$k == 2, ]
  both <- !is.na(a$allele) & !is.na(a$other) & !is.na(b$allele) & !is.na(b$other)
  same <- both & b$allele == a$allele & b$other == a$other
  swap <- both & b$allele == a$other & b$other == a$allele
  bad <- both & !same & !swap & is.finite(b$beta)
  b$beta[swap] <- -b$beta[swap]
  if (any(bad)) {
    message("Alleles of ", sum(bad), " variant(s) differ between the two GWAS (strand ",
            "flips or different variants), so ", names[2], "'s estimate is left out: ",
            paste(reps$snp[bad], collapse = ", "), ".")
    b$beta[bad] <- NA
  }
  unsure <- !both & is.finite(a$beta) & is.finite(b$beta) &
    pmax(a$maf, b$maf, na.rm = TRUE) > 0.4
  unsure[is.na(unsure)] <- FALSE
  if (any(unsure)) {
    message("Without allele columns in both GWAS, effects are compared per each GWAS's ",
            "minor allele (or as reported); these variants have a minor allele frequency ",
            "above 40%, so their minor alleles may differ: ",
            paste(reps$snp[unsure], collapse = ", "), ".")
  }
  est[est$k == 2, ] <- b
  est
}

.short_id <- function(x) {
  sub("^((chr)?[0-9XYMT]+[:_][0-9]+)[:_][^:_]+[:_][^:_]+$", "\\1", x, ignore.case = TRUE)
}
