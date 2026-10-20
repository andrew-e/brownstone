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
#' * **Effects that disagree with their own allele frequencies are flagged.**
#'   For binary traits whose summary statistics report the allele frequency
#'   among cases and among controls, a lead variant whose reported effect
#'   disagrees with the crude odds ratio those frequencies imply (by an order
#'   of magnitude, or by six standard errors) is marked with a black ring; see
#'   [find_artifacts()].
#' * **The detection boundary** is the smallest effect this study could have
#'   found at each frequency, drawn as a curve beneath the cloud: a true effect
#'   below it is invisible to this study, whatever the biology. The typical
#'   standard error grows as `c / sqrt(2f(1 - f))` with frequency *f*, for
#'   quantitative traits and log odds ratios alike, and its constant *c* (the
#'   sample size, case/control balance and test calibration of the study as
#'   actually run, not an idealised one) is estimated from the variants
#'   themselves: the median of `SE x sqrt(2f(1 - f))`. The boundary is then
#'   *(z*(1-*α*/2)* + z*(power)*) x typical SE*, with *α* the `threshold`.
#'   `detection = FALSE` removes it; `power` sets the level, or levels.
#'
#' Variants without an allele frequency, or with a frequency of 0, are left out.
#' Effects are not corrected for winner's curse, which inflates rare variants
#' the most.
#'
#' @inheritParams volcano_plot
#' @param effects `"magnitude"` (the default): |beta|, or the odds ratio of the
#'   risk allele. `"signed"`: beta, oriented to the minor allele where an allele
#'   frequency column says which allele that is.
#' @param detection Draw the detection boundary (see above).
#' @param power Power of the detection boundary: a number between 0 and 1, or
#'   several, one curve each.
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
                          point_size = 1, detection = TRUE, power = 0.8) {
  build <- match.arg(build)
  effects <- match.arg(effects)
  effect_scale <- match.arg(effect_scale)
  if (any(power <= 0 | power >= 1, na.rm = TRUE)) {
    stop("`power` is a power, not a percentage: give it between 0 and 1 (e.g. 0.8).",
         call. = FALSE)
  }
  signed <- effects == "signed"
  label_col <- if (is.character(labels) && length(labels) == 1L &&
                    !labels %in% c("nearest_gene", "snp")) labels
  d <- .standardise(data, ..., keep = label_col, ids = FALSE, ccafs = TRUE)
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
  # The typical standard error's constant, from every variant that has an SE,
  # before any thinning (see .detection_boundary()).
  c_se <- if (detection) .se_unit(d)
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
  # Lead variants whose effects disagree with their own case/control allele
  # frequencies (see find_artifacts()), marked with a ring among the hits.
  d$flag <- seq_len(nrow(d)) %in% .artifact_flags(d, threshold, window, build)

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

  # The detection boundary, over the drawn frequency range.
  bnd <- if (detection) .detection_boundary(c_se, lo, threshold, power, signed)

  ggplot2::ggplot() +
    ggplot2::annotate("rect", xmin = bands$xmin, xmax = bands$xmax, ymin = -Inf, ymax = Inf,
                      fill = .ink$band, colour = NA) +
    ggplot2::geom_hline(yintercept = 0, colour = .ink$axis, linewidth = 0.4) +
    ggplot2::geom_point(data = noise, ggplot2::aes(x = .data$x, y = .data$y),
                        colour = .ink$noise, size = 0.45 * pts, stroke = 0) +
    .boundary_layers(bnd, base_size) +
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

# The constant of the typical standard error: SE(beta) ~ c / sqrt(2 f (1 - f))
# with frequency f, for quantitative traits and log odds ratios alike. c is
# estimated as the median of SE x sqrt(2 f (1 - f)) over the variants that have
# an SE, so it absorbs the sample size, the case/control balance and the
# calibration of the test as actually run, rather than an idealised study.
.se_unit <- function(d) {
  ok <- is.finite(d$se) & is.finite(d$maf) & d$maf > 0
  if (!any(ok)) return(NA_real_)
  stats::median(d$se[ok] * sqrt(2 * d$maf[ok] * (1 - d$maf[ok])))
}

# The detection boundary over the frequency range [lo, 0.5]: the smallest
# effect this study could detect at each frequency, at each of the given
# powers, (z(1 - alpha/2) + z(power)) x c / sqrt(2 f (1 - f)), with alpha the
# significance threshold. Mirrored about zero for signed effects. One curve
# per power (and side), distinguished by `curve`.
.detection_boundary <- function(c_se, lo, threshold, power, signed) {
  if (!is.finite(c_se)) return(NULL)
  f <- 10^seq(log10(lo), log10(0.5), length.out = 200)
  se <- c_se / sqrt(2 * f * (1 - f))   # the typical SE at each frequency
  sides <- if (signed) c(1, -1) else 1
  curves <- unlist(
    lapply(sort(power), function(p) {
      z <- stats::qnorm(1 - threshold / 2) + stats::qnorm(p)
      lapply(sides, function(s) {
        data.frame(curve = paste(p, s, sep = "-"), power = p, side = s,
                   x = log10(f), y = s * z * se, stringsAsFactors = FALSE)
      })
    }),
    recursive = FALSE
  )
  do.call(rbind, curves)
}

# The boundary's layers: the curve(s) over the cloud, and one label per power
# at the common (right-hand) end of its curve, where it sits lowest.
.boundary_layers <- function(bnd, base_size) {
  if (is.null(bnd)) return(NULL)
  up <- bnd[bnd$side == 1, ]
  lab <- up[!duplicated(up$power, fromLast = TRUE), ]
  list(
    ggplot2::geom_line(
      data = bnd, ggplot2::aes(x = .data$x, y = .data$y, group = .data$curve),
      colour = .ink$secondary, linewidth = 0.35
    ),
    ggplot2::annotate(
      "text", x = lab$x, y = lab$y, hjust = 1, vjust = -0.5,
      label = sprintf("%g%% power", 100 * lab$power),
      size = base_size * .text_size$note / ggplot2::.pt, colour = .ink$muted
    )
  )
}
