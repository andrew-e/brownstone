#' Volcano plot: effect size against significance
#'
#' @description
#' Plots each variant's effect size against its -log10(P), in the style of
#' [brownstone()]: a grey cloud for the bulk of the genome, colour for the hits.
#'
#' * **Effects are signed**, oriented to the minor allele where an allele
#'   frequency is available (see `effect_allele` in [standardise_sumstats()]),
#'   so that the sign means the same thing at every locus. Odds ratios are drawn
#'   on the log scale and labelled as odds ratios.
#' * **Hits are coloured by minor allele frequency**, as in [brownstone()]:
#'   variants past the significance line are solid, peaks that only reach the
#'   `suggestive` line are hollow rings, and everything else is a grey cloud,
#'   thinned where points would be drawn on top of each other.
#' * **The effect axis reacts to the data.** Rare variants have large standard
#'   errors, so even null ones can have large estimates. The axis is linear
#'   unless a few estimates dwarf the lead effects; then, as on the effect track
#'   of [brownstone()], it breaks and puts the outliers on a log scale
#'   (`effect_scale`).
#' * **Loci are labelled** with their nearest gene (or `labels`), one label per
#'   locus, most significant first.
#' * **Effects that disagree with their own allele frequencies are flagged.**
#'   For binary traits whose summary statistics report the allele frequency
#'   among cases and among controls, a lead variant whose reported effect
#'   disagrees with the crude odds ratio those frequencies imply (by an order
#'   of magnitude, or by six standard errors) is marked with a black ring; see
#'   [find_artifacts()].
#'
#' @inheritParams brownstone
#' @param labels How to label loci, as in [brownstone()]: `"nearest_gene"`,
#'   `"snp"`, a column name, a data frame of labels, or `NULL` for none.
#' @param max_labels Label at most this many loci (the most significant).
#' @param effect_scale How the effect axis reacts to the data: `"auto"` (the
#'   default; linear, breaking for outlying estimates), `"sqrt"` or `"linear"`.
#'   See [brownstone()].
#'
#' @return A ggplot.
#' @export
#' @examples
#' gwas <- simulate_gwas(n_variants = 2e5)
#' volcano_plot(gwas)
volcano_plot <- function(data, ..., build = c("GRCh38", "GRCh37"), threshold = 5e-8,
                         suggestive = NULL, labels = "nearest_gene", max_labels = 20,
                         window = 5e5, effect_scale = c("auto", "sqrt", "linear"),
                         maf_breaks = c(0.001, 0.01, 0.05), colours = NULL,
                         y_compress = "auto", base_size = 10, point_size = 1) {
  build <- match.arg(build)
  effect_scale <- match.arg(effect_scale)
  label_col <- if (is.character(labels) && length(labels) == 1L &&
                    !labels %in% c("nearest_gene", "snp")) labels
  d <- .standardise(data, ..., keep = label_col, ids = FALSE, ccafs = TRUE)
  if (!isTRUE(attr(d, "has_effect"))) {
    stop("volcano_plot() needs effect size and standard error columns; ",
         "name them with `beta =` (or `or =`) and `se =`.", call. = FALSE)
  }
  sig <- -log10(threshold)
  if (!is.null(suggestive) && suggestive <= threshold) {
    stop("`suggestive` must be a larger P-value than `threshold`.", call. = FALSE)
  }
  sugg <- if (!is.null(suggestive)) -log10(suggestive)
  floor_logp <- sugg %||% sig
  has_maf <- isTRUE(attr(d, "has_maf"))
  if (!has_maf) message("No allele frequency column found; hits are drawn in one colour.")
  log_or <- identical(attr(d, "effect_scale"), "log_or")
  allele <- attr(d, "effect_allele")

  d <- d[is.finite(d$beta), ]   # a variant with no estimate has no place on the x axis
  d$x <- d$beta
  max_logp <- max(d$logp)
  # A table of labels can name any variant, so it is matched before thinning.
  if (is.data.frame(labels)) labels <- .label_table(labels, d)

  # The effect axis is set by the lead effects, as on brownstone()'s track.
  cl <- .clump(d, floor_logp, window)
  lead_abs <- abs(d$beta[cl$lead[d$logp[cl$lead] >= sig]])
  top <- max(abs(d$beta))
  ax <- .effect_axis(lead_abs, top, effect_scale, max_estimate = top)

  # Thin the grey cloud on the axis as drawn, then work on far fewer rows.
  # Every hit survives (.thin() only drops points below `floor_logp`), so the
  # lead variants are carried through as a flag.
  keep <- .thin(ax$trans$transform(d$x), d$logp, below = floor_logp)
  d$lead <- FALSE
  d$lead[cl$lead] <- TRUE
  d <- d[keep, ]
  leads <- which(d$lead)[order(-d$logp[d$lead])]   # most significant first
  # Lead variants whose effects disagree with their own case/control allele
  # frequencies (see find_artifacts()), marked with a ring among the hits.
  d$flag <- seq_len(nrow(d)) %in% .artifact_flags(d, threshold, window, build)

  cp <- .class_palette(d, has_maf, maf_breaks, colours, floor_logp)
  d$class <- cp$class
  cs <- .colour_spec(d, has_maf, "binned", cp$pal, NULL, floor_logp)
  d$col <- cs$values
  d$stack <- cs$stack

  comp <- .y_compressor(max_logp, sig, y_compress)
  d$y <- comp$f(d$logp)
  # The cloud stays grey. Tinting it by frequency class shows each class's fan,
  # but tints calm enough not to compete with the hits can't be told apart with
  # colour blindness (deutan delta-E 3.9 between the aqua and orange tints).
  noise <- d[d$logp < floor_logp, ]
  hits <- d[d$logp >= floor_logp, ]
  hits$tier <- factor(ifelse(hits$logp >= sig, "gws", "suggestive"),
                      levels = c("gws", "suggestive"))
  ymax <- max(comp$f(max_logp) * 1.1, sig * 1.6)
  pts <- point_size * base_size / 10

  lab <- .locus_labels(labels, label_col, d, leads, sig, build, max_labels)
  lab$y <- comp$f(lab$logp)
  near <- hits[.thin(ax$trans$transform(hits$x), hits$y, below = Inf, nx = 600,
                     ny_per_unit = 1), ]

  p <- ggplot2::ggplot()
  if (!is.null(comp$at)) p <- p + .axis_break(comp$at, ymax)
  p +
    ggplot2::geom_vline(xintercept = 0, colour = .ink$axis, linewidth = 0.4) +
    ggplot2::geom_point(data = noise, ggplot2::aes(x = .data$x, y = .data$y),
                        colour = .ink$noise, size = 0.45 * pts, stroke = 0) +
    .effect_break(ax, top, signed = TRUE, axis = "x") +
    .p_lines(sig, sugg) +
    .hit_layers(hits, cs, legend = has_maf, sig, sugg, pts) +
    .repel_labels(lab, near, base_size, nudge_y = 0.04 * ymax,
                  ylim = c(sig + 0.03 * ymax, NA)) +
    ggplot2::scale_x_continuous(
      transform = ax$trans, limits = c(-top, top),
      breaks = .track_breaks(ax, top, signed = TRUE, log_or = log_or),
      labels = function(b) .as_plotmath(.effect_labels(b, log_or)),
      expand = ggplot2::expansion(mult = 0.03)
    ) +
    .p_axis(comp, max_logp, ymax, sig, sugg) +
    ggplot2::labs(x = .effect_title(signed = TRUE, log_or, allele),
                  y = expression(-log[10](italic(P)))) +
    theme_brownstone(base_size = base_size) + .legend_on_top(base_size) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(colour = .ink$grid, linewidth = 0.3),
      axis.text.y.right = ggplot2::element_text(colour = .ink$muted, size = base_size * .text_size$note,
                                                hjust = 0)
    )
}
