#' Effect size against minor allele frequency
#'
#' @description
#' Plots each variant's effect size against its minor allele frequency, the
#' view that shows how rare variants with large effects and common variants
#' with small ones fit together in one GWAS.
#'
#' * **Frequency is on a log scale**, which spreads out the rare end where most
#'   variants are. Alternating bands mark the frequency classes of [brownstone()]
#'   (`maf_breaks`).
#' * **Significance is shown by style, not colour** (frequency is already on the
#'   x axis): genome-wide significant variants are solid, peaks that only reach
#'   the `suggestive` line are hollow rings, and everything else is a grey cloud,
#'   thinned where points would be drawn on top of each other.
#' * **Effects are |beta| by default** (or the odds ratio of the risk allele):
#'   the size of an effect, whichever allele the GWAS reported. With
#'   `effects = "signed"`, beta is oriented to the minor allele. As on the effect
#'   track of [brownstone()], the axis breaks for a few outlying estimates
#'   (`effect_scale`).
#' * **Genome-wide significant loci are labelled** with their nearest gene (or
#'   `labels`), at their lead variant.
#'
#' Variants without an allele frequency, or with a frequency of 0, are left out.
#' Effects are not corrected for winner's curse, which inflates rare variants
#' the most.
#'
#' @inheritParams volcano_plot
#' @param effects `"magnitude"` (the default): |beta|, or the odds ratio of the
#'   risk allele. `"signed"`: beta, oriented to the minor allele where an allele
#'   frequency column says which allele that is.
#' @param maf_breaks Boundaries of the minor allele frequency classes, marked by
#'   alternating bands and axis ticks.
#'
#' @return A ggplot.
#' @export
#' @examples
#' gwas <- simulate_gwas(n_variants = 2e5)
#' effect_vs_maf(gwas)
effect_vs_maf <- function(data, ..., build = c("GRCh38", "GRCh37"),
                          effects = c("magnitude", "signed"), threshold = 5e-8,
                          suggestive = NULL, labels = "nearest_gene", max_labels = 10,
                          window = 5e5, effect_scale = c("auto", "sqrt", "linear"),
                          maf_breaks = c(0.001, 0.01, 0.05), base_size = 10,
                          point_size = 1) {
  build <- match.arg(build)
  effects <- match.arg(effects)
  effect_scale <- match.arg(effect_scale)
  signed <- effects == "signed"
  label_col <- if (is.character(labels) && length(labels) == 1L &&
                   !labels %in% c("nearest_gene", "snp")) labels
  d <- .standardise(data, ..., keep = label_col, ids = FALSE)
  if (!isTRUE(attr(d, "has_effect"))) {
    stop("effect_vs_maf() needs effect size and standard error columns; ",
         "name them with `beta =` (or `or =`) and `se =`.", call. = FALSE)
  }
  if (!isTRUE(attr(d, "has_maf"))) {
    stop("effect_vs_maf() needs an allele frequency column; ",
         "name it with `maf =` or `eaf =`.", call. = FALSE)
  }
  sig <- -log10(threshold)
  if (!is.null(suggestive) && suggestive <= threshold) {
    stop("`suggestive` must be a larger P-value than `threshold`.", call. = FALSE)
  }
  sugg <- if (!is.null(suggestive)) -log10(suggestive)
  floor_logp <- sugg %||% sig
  log_or <- identical(attr(d, "effect_scale"), "log_or")
  allele <- attr(d, "effect_allele")

  d <- d[is.finite(d$beta) & is.finite(d$maf) & d$maf > 0, ]
  d$y <- if (signed) d$beta else abs(d$beta)
  # Rows are tracked by index while labels are matched (a table of labels can
  # name any variant, so it is matched before thinning).
  d$x <- seq_len(nrow(d))
  if (is.data.frame(labels)) labels <- .label_table(labels, d)
  cl <- .clump(d, floor_logp, window)
  gws_leads <- cl$lead[d$logp[cl$lead] >= sig]
  lab <- .locus_labels(labels, label_col, d, gws_leads, sig, build, max_labels)
  lab$y <- d$y[lab$x]
  lab$x <- log10(d$maf[lab$x])

  # The effect axis is set by the lead effects, as on brownstone()'s track.
  top <- max(abs(d$y))
  ax <- .effect_axis(abs(d$beta[gws_leads]), top, effect_scale, max_estimate = top)
  lim <- if (signed) c(-top, top) else c(0, top)

  # Frequency goes on the x axis as log10(MAF), labelled in percent. Thin the
  # grey cloud on the axes as drawn.
  d$x <- log10(d$maf)
  low <- which(d$logp < floor_logp)
  yt <- ax$trans$transform(d$y[low])
  ys <- (yt - min(yt)) / max(diff(range(yt)), 1e-12) * 1500
  keep <- c(low[.thin(d$x[low], ys, below = Inf, nx = 3000, ny_per_unit = 1)],
            which(d$logp >= floor_logp))
  d <- d[sort(keep), ]

  noise <- d[d$logp < floor_logp, ]
  hits <- d[d$logp >= floor_logp, ]
  hits$tier <- factor(ifelse(hits$logp >= sig, "gws", "suggestive"),
                      levels = c("gws", "suggestive"))
  # One colour for the hits: frequency is the x axis.
  hits$class <- factor(rep("all", nrow(hits)))
  cs <- .colour_spec(hits, FALSE, "binned", c(all = maf_colours(1)), NULL, floor_logp)
  hits$col <- cs$values
  hits$stack <- cs$stack
  near <- hits[.thin(hits$x, ax$trans$transform(hits$y), below = Inf, nx = 600,
                     ny_per_unit = 30), ]
  pts <- point_size * base_size / 10

  # Frequency classes as alternating bands, rarest at the left.
  lo <- 10^floor(min(d$x))
  edges <- log10(sort(unique(c(lo, maf_breaks[maf_breaks > lo & maf_breaks < 0.5], 0.5))))
  bands <- data.frame(xmin = edges[-length(edges)], xmax = edges[-1])
  bands <- bands[seq_len(nrow(bands)) %% 2 == 0, ]
  breaks <- sort(unique(c(10^seq(log10(lo), -1), maf_breaks, 0.5)))
  breaks <- breaks[breaks >= lo & breaks <= 0.5]

  ggplot2::ggplot() +
    ggplot2::annotate("rect", xmin = bands$xmin, xmax = bands$xmax, ymin = -Inf, ymax = Inf,
                      fill = .ink$band, colour = NA) +
    ggplot2::geom_hline(yintercept = 0, colour = .ink$axis, linewidth = 0.4) +
    ggplot2::geom_point(data = noise, ggplot2::aes(x = .data$x, y = .data$y),
                        colour = .ink$noise, size = 0.45 * pts, stroke = 0) +
    .effect_break(ax, top, signed = signed, axis = "y") +
    .hit_layers(hits, cs, legend = FALSE, sig, sugg, pts) +
    .repel_labels(lab, near, base_size) +
    ggplot2::scale_x_continuous(
      limits = log10(c(lo, 0.5)), breaks = log10(breaks), expand = ggplot2::expansion(mult = 0.02),
      labels = paste0(format(breaks * 100, drop0trailing = TRUE, trim = TRUE,
                             scientific = FALSE), "%")
    ) +
    ggplot2::scale_y_continuous(
      transform = ax$trans, limits = lim, oob = scales::oob_squish,
      breaks = .track_breaks(ax, top, signed = signed, log_or = log_or),
      labels = function(b) .as_plotmath(.effect_labels(b, log_or)),
      expand = ggplot2::expansion(mult = c(if (signed) 0.04 else 0.01, 0.04))
    ) +
    ggplot2::labs(x = "Minor allele frequency", y = .effect_title(signed, log_or, allele)) +
    theme_brownstone(base_size = base_size) + .legend_on_top(base_size) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}
