#' Miami plot: two Manhattan plots, one above the other
#'
#' @description
#' Compares two GWAS (for example two sexes, two ancestries or two related
#' traits) on one genome axis. The first GWAS's Manhattan rises above the
#' chromosome axis and the second's hangs below it; beneath them, an effect
#' track does the same for |beta|: the first GWAS's effects rise from a zero
#' line and the second's hang from it. Everything else is as in
#' [manhattan()]: hits are coloured by minor allele frequency (one palette for
#' both GWAS), loci are labelled with their nearest gene, and both axes break
#' to make room for very strong signals.
#'
#' * **One scale for both.** The two halves share the -log10(P) axis (and its
#'   break) and the |beta| axis, so a taller tower or a longer effect means the
#'   same thing above and below the line. When one GWAS is much stronger, the
#'   other looks flat: that is the comparison.
#' * **Magnitudes only.** The effect track shows |beta| (or the odds ratio of
#'   the risk allele), so the two GWAS can be compared without aligning their
#'   alleles. For signed effects side by side, see `forest_plot(compare = )`.
#' * **Different file formats.** Column names given in `...` apply to both
#'   GWAS. When the two files name their columns differently, standardise each
#'   first with [standardise_sumstats()] and pass the results.
#'
#' @inheritParams manhattan
#' @param top,bottom GWAS summary statistics, one row per variant, for the
#'   upper and lower halves. Columns are detected automatically; see
#'   [standardise_sumstats()].
#' @param ... Column names and options passed to [standardise_sumstats()] for
#'   both GWAS.
#' @param names Names of the two GWAS, printed beside each half: top left of
#'   the upper half, bottom left of the lower half. By default, the expressions
#'   passed as `top` and `bottom`.
#' @param labels How to label loci, as in [manhattan()]. A data frame of labels
#'   labels each listed variant in whichever GWAS it is found.
#' @param max_labels Label at most this many loci in each half (the most
#'   significant).
#' @param effects `"magnitude"` (the default) draws the |beta| track; `"none"`
#'   (or `FALSE`) leaves it out. Both GWAS need effect and standard error
#'   columns, on the same scale (both beta, or both odds ratios).
#' @param track_height Fraction of the figure height given to the effect track
#'   (both halves together).
#'
#' @return A patchwork object. Add a title with
#'   `+ patchwork::plot_annotation(title = "...")`.
#' @export
#' @examples
#' men <- simulate_gwas(n_variants = 2e5, seed = 1)
#' women <- simulate_gwas(n_variants = 2e5, seed = 2)
#' miami_plot(men, women)
miami_plot <- function(top, bottom, ..., names = NULL, build = c("GRCh38", "GRCh37"),
                       threshold = 5e-8, suggestive = NULL,
                       labels = "nearest_gene", label_suggestive = TRUE,
                       track_suggestive = TRUE, max_labels = 25, window = 5e5,
                       effects = c("magnitude", "none"),
                       effect_scale = c("auto", "sqrt", "linear"),
                       allelic_series = TRUE,
                       maf_scale = c("binned", "continuous"),
                       maf_breaks = c(0.001, 0.01, 0.05), colours = NULL,
                       y_compress = "auto", track_height = 0.3,
                       track_variants = c("all", "leads"),
                       label_angle = 0, base_size = 10, point_size = 1) {
  names <- names %||% c(deparse1(substitute(top)), deparse1(substitute(bottom)))
  if (length(names) != 2L) stop("`names` must name the two GWAS.", call. = FALSE)
  names <- as.character(names)
  build <- match.arg(build)
  effect_scale <- match.arg(effect_scale)
  track_variants <- match.arg(track_variants)
  maf_scale <- match.arg(maf_scale)
  if (isTRUE(effects)) effects <- "magnitude"
  if (is.null(effects) || isFALSE(effects)) effects <- "none"
  effects <- match.arg(effects)
  label_col <- if (is.character(labels) && length(labels) == 1L &&
                   !labels %in% c("nearest_gene", "snp")) labels
  sig <- -log10(threshold)
  if (!is.null(suggestive) && suggestive <= threshold) {
    stop("`suggestive` must be a larger P-value than `threshold`.", call. = FALSE)
  }
  sugg <- if (!is.null(suggestive)) -log10(suggestive)
  floor_logp <- sugg %||% sig

  d <- list(.standardise(top, ..., keep = label_col, ids = FALSE),
            .standardise(bottom, ..., keep = label_col, ids = FALSE))
  # One genome axis for both, so that a locus sits at the same x in each half
  # even when one GWAS lacks a chromosome.
  layout <- .genome_layout(data.frame(chr = c(d[[1]]$chr, d[[2]]$chr),
                                      pos = c(d[[1]]$pos, d[[2]]$pos)))
  max_logp <- max(d[[1]]$logp, d[[2]]$logp)
  if (is.data.frame(labels)) {
    found <- !is.na(.match_variants(labels, d[[1]])) | !is.na(.match_variants(labels, d[[2]]))
    if (!all(found)) {
      warning(sum(!found), " labelled variant(s) not found in either GWAS.", call. = FALSE)
    }
  }
  half_labels <- vector("list", 2L)
  for (k in 1:2) {
    d[[k]]$x <- .genome_x(d[[k]]$chr, d[[k]]$pos, layout)
    # A table of labels can name any variant, so it is matched before thinning.
    half_labels[[k]] <- if (is.data.frame(labels)) .label_table(labels, d[[k]], warn = FALSE) else labels
    d[[k]] <- d[[k]][.thin(d[[k]]$x, d[[k]]$logp, below = floor_logp), ]
  }

  # One palette for both halves: classes, an "Unknown" class and the range of a
  # continuous scale are set by the two GWAS together.
  has_maf <- isTRUE(attr(d[[1]], "has_maf")) || isTRUE(attr(d[[2]], "has_maf"))
  n1 <- nrow(d[[1]])
  both <- data.frame(maf = c(d[[1]]$maf, d[[2]]$maf), logp = c(d[[1]]$logp, d[[2]]$logp))
  cset <- .colour_setup(both, has_maf, maf_scale, maf_breaks, colours, floor_logp)
  cs <- cset$cs
  loc <- vector("list", 2L)
  for (k in 1:2) {
    rows <- if (k == 1) seq_len(n1) else n1 + seq_len(nrow(d[[2]]))
    for (v in c("class", "col", "stack")) d[[k]][[v]] <- cset[[v]][rows]
    loc[[k]] <- .locus_setup(d[[k]], floor_logp, sig, window, half_labels[[k]], label_col,
                             label_suggestive, build, max_labels)
    d[[k]] <- loc[[k]]$d
  }
  comp <- .y_compressor(max_logp, sig, y_compress)

  has_effects <- effects != "none"
  if (has_effects) {
    eff_ok <- vapply(d, function(x) isTRUE(attr(x, "has_effect")), logical(1))
    scales <- vapply(d, function(x) attr(x, "effect_scale") %||% "beta", "")
    on_track <- vapply(loc, function(l) any(l$strong_locus | (track_suggestive & !is.null(sugg))),
                       logical(1))
    if (!all(eff_ok)) {
      message("No effect size and standard error columns found in ",
              paste(names[!eff_ok], collapse = " or "), "; drawing without the effect track.")
      has_effects <- FALSE
    } else if (scales[1] != scales[2]) {
      message("One GWAS reports odds ratios and the other betas, which can't share an axis; ",
              "drawing without the effect track.")
      has_effects <- FALSE
    } else if (!any(on_track)) {
      has_effects <- FALSE
    }
  }

  # The legend goes with the upper half, unless only the lower half has hits.
  hits_in <- vapply(d, function(x) any(x$logp >= floor_logp), logical(1))
  legend_top <- hits_in[1] || !hits_in[2]
  upper <- .main_panel(d[[1]], max_logp, layout, comp, sig, sugg, loc[[1]]$lab, cs, has_maf,
                       label_angle, base_size, point_size, x_axis = TRUE, floor_gap = 0.03,
                       legend = legend_top, top_title = names[1])
  lower <- .main_panel(d[[2]], max_logp, layout, comp, sig, sugg, loc[[2]]$lab, cs, has_maf,
                       label_angle, base_size, point_size, x_axis = FALSE, floor_gap = 0.03,
                       down = TRUE, legend = !legend_top)
  # The chromosome names run between the two halves.
  upper <- upper + .half_names(NULL, NULL, base_size) + ggplot2::labs(x = NULL) +
    ggplot2::theme(plot.margin = ggplot2::margin(6, 10, 0, 6))
  lower <- lower + .half_names(NULL, names[2], base_size) +
    ggplot2::theme(plot.margin = ggplot2::margin(2, 10, if (has_effects) 0 else 4, 6))
  # The legend is collected above the figure, so that nothing comes between
  # the name of the upper half and its panel; the two halves share one
  # -log10(P) axis title.
  if (!has_effects) {
    return(patchwork::wrap_plots(upper, lower, ncol = 1, guides = "collect",
                                 axis_titles = "collect_y") +
             patchwork::plot_annotation(theme = theme_maffect(base_size = base_size)))
  }

  tr <- lapply(1:2, function(k) {
    .track_setup(d[[k]], loc[[k]]$cl, sig, floor_logp, window, allelic_series && has_maf,
                 track_suggestive, track_variants)
  })
  eff <- rbind(cbind(tr[[1]]$eff, side = rep(1, nrow(tr[[1]]$eff))),
               cbind(tr[[2]]$eff, side = rep(-1, nrow(tr[[2]]$eff))))
  cloud <- if (track_variants == "all") {
    rbind(cbind(tr[[1]]$cloud, side = rep(1, nrow(tr[[1]]$cloud))),
          cbind(tr[[2]]$cloud, side = rep(-1, nrow(tr[[2]]$cloud))))
  }
  track <- .effect_panel(eff, cloud, layout, cs, mode = "mirror", scale = effect_scale,
                         log_or = identical(attr(d[[1]], "effect_scale"), "log_or"),
                         allele = NULL, down = FALSE, base_size = base_size,
                         point_size = point_size) +
    .half_names(names[1], names[2], base_size) + ggplot2::labs(x = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(10, 10, 4, 6))
  # A bold rule separates the two mirrored plots: -log10(P) above, |beta|
  # below.
  rule <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, colour = .ink$secondary, linewidth = .ground_width) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = ggplot2::expansion(0)) +
    ggplot2::theme_void() + ggplot2::theme(plot.margin = ggplot2::margin(4, 10, 4, 6))
  h <- (1 - track_height) / 2
  patchwork::wrap_plots(upper, lower, rule, track, ncol = 1,
                        heights = c(h, h, 0.01, track_height),
                        guides = "collect", axis_titles = "collect_y") +
    patchwork::plot_annotation(theme = theme_maffect(base_size = base_size))
}

# Each half of a Miami panel is named in the margin, top left for the upper
# half and bottom left for the lower one (a caption): inside the panel a name
# would collide with towers and their labels. The upper -log10(P) panel is
# named by its top x-axis title instead (see .x_scale()), which, unlike a
# title, sits below the legend, right against its panel.
.half_names <- function(upper, lower, base_size) {
  style <- ggplot2::element_text(colour = .ink$primary, face = "bold", hjust = 0,
                                 size = base_size * .text_size$legend,
                                 margin = ggplot2::margin(2, 0, 2, 0))
  list(ggplot2::labs(title = upper, caption = lower),
       ggplot2::theme(plot.title = style, plot.caption = style, axis.title.x.top = style,
                      plot.title.position = "panel", plot.caption.position = "panel"))
}
