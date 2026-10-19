#' QQ plot of GWAS P-values, by allele frequency, with genomic inflation
#'
#' @description
#' Plots each variant's observed -log10(P) against the value expected if no
#' variant were associated. By default there is one curve per minor allele
#' frequency class, in the same colours as [brownstone()], each compared with
#' its own expected distribution: test statistics for rare variants are often
#' calibrated differently from common ones, and a single curve hides that.
#' The genomic inflation factor, \eqn{\lambda_{GC}}, is printed on the plot,
#' for all variants and for each class.
#'
#' * **Genomic inflation** is the median of the 1-df chi-square statistics
#'   implied by the P-values, divided by its expected value (0.455). It is
#'   computed from every variant, before any thinning for plotting. It assumes
#'   1-df tests: gene-based or 2-df tests break it. Under polygenicity the
#'   common-variant \eqn{\lambda} grows with sample size even without
#'   confounding (the LD score regression intercept tells the two apart), and
#'   P-values already corrected by genomic control give \eqn{\lambda \approx 1}
#'   by construction.
#' * **Rare variants** are often conservative (\eqn{\lambda < 1}, a curve below
#'   the diagonal) with saddlepoint (SPA) or Firth tests of binary traits, and
#'   with low minor allele counts. \eqn{\lambda > 1} in a rare class of an
#'   unbalanced binary trait tested without SPA suggests miscalibration.
#' * **The grey bands** are pointwise 95% intervals for each class's ordered
#'   P-values under the null (from the Beta distribution of uniform order
#'   statistics), one per class: smaller classes have wider bands that end
#'   earlier.
#' * **Other strata.** `by` can also name a column, for example `"chr"`, or a
#'   binned imputation quality score kept with `keep =`. Up to four groups are
#'   drawn as coloured curves on one plot; more get a panel each, with their
#'   \eqn{\lambda} in the panel title, because more than four colours cannot
#'   all stay distinguishable with colour blindness.
#'
#' Below the significance line, points that would be drawn on top of each
#' other are thinned (as in [brownstone()]), so millions of variants plot
#' quickly.
#'
#' @inheritParams brownstone
#' @param by `"maf"` (the default) draws one curve per minor allele frequency
#'   class (see `maf_breaks`); `"none"` draws one curve for all variants. Without
#'   an allele frequency column there is one curve. Any other value names a
#'   column of the standardised data (see [standardise_sumstats()]): `"chr"`,
#'   or a column of `data` kept with `keep = "<column>"`. There is one curve (or
#'   panel) per value, and variants with no value form an "Unknown" group. A
#'   continuous column should be binned first: at most 30 groups are drawn.
#' @param colours One colour per frequency class (or group), most common first
#'   (default [maf_colours()]); for `by = "none"`, the colour of the one curve.
#' @param lambda Print \eqn{\lambda_{GC}} on the plot: for all variants in the
#'   corner of the panel (or the subtitle, with a panel per group) and, for
#'   classes of at least 1,000 variants, in the legend (or panel titles).
#' @param ci Draw the pointwise 95% bands.
#' @param threshold Genome-wide significance threshold, drawn as a line.
#'
#' @return A ggplot.
#' @export
#' @examples
#' gwas <- simulate_gwas(n_variants = 2e5)
#' qq_plot(gwas)
qq_plot <- function(data, ..., by = c("maf", "none"), maf_breaks = c(0.001, 0.01, 0.05),
                    colours = NULL, lambda = TRUE, ci = TRUE, threshold = 5e-8,
                    y_compress = "auto", base_size = 10, point_size = 1) {
  by <- by[1]   # "maf" by default
  if (!is.character(by) || is.na(by)) {
    stop("`by` must be \"maf\", \"none\" or a column name.", call. = FALSE)
  }
  d <- .standardise(data, ..., ids = FALSE)
  if (!by %in% c("maf", "none", names(d))) {
    stop(sprintf(paste0("`by` must be \"maf\", \"none\" or a column of the standardised ",
                        "data (%s). To use another column of `data`, keep it too: ",
                        "`keep = \"%s\"`."), paste(names(d), collapse = ", "), by),
         call. = FALSE)
  }
  sig <- -log10(threshold)
  logp <- d$logp
  if (by == "maf" && !isTRUE(attr(d, "has_maf"))) {
    message("No allele frequency column found; drawing one curve for all variants.")
    by <- "none"
  }

  # A QQ plot assumes every tested variant is there: about half the P-values
  # should be above 0.5 (still about 42% with lambda = 1.5).
  above_half <- mean(logp < log10(2))
  if (above_half < 0.25) {
    warning(sprintf(paste0("Only %.0f%% of P-values are above 0.5 (about 50%% expected): ",
                           "the data look filtered. The QQ plot assumes every tested variant ",
                           "is included; genomic inflation is not shown."), 100 * above_half),
            call. = FALSE)
    lambda <- FALSE
  }
  at_one <- mean(logp == 0)
  if (at_one >= 0.01) {
    warning(sprintf(paste0("%.1f%% of P-values are exactly 1 (failed or degenerate tests?); ",
                           "they deflate genomic inflation."), 100 * at_one), call. = FALSE)
  }

  # Strata: frequency classes or the values of a column, with variants of
  # unknown frequency (or value) as their own class rather than dropped. Up to
  # four strata are told apart by colour (the frequency-class palette, which
  # stays distinguishable with colour blindness); more get a panel each.
  facet <- FALSE
  if (by == "none") {
    cls <- factor(rep("all", length(logp)))
    pal <- c(all = (colours %||% maf_colours(1))[1])
  } else {
    v <- if (by == "maf") NULL else d[[by]]
    cls <- if (by == "maf") {
      .maf_classes(d$maf, maf_breaks)
    } else if (by == "chr") {
      factor(v, intersect(.chr_levels, unique(v)))
    } else if (is.factor(v)) {
      droplevels(v)
    } else {
      factor(v)
    }
    n <- nlevels(cls)
    if (n > 30) {
      stop(sprintf("`by = \"%s\"` gives %d groups, too many to draw; bin a continuous column first.",
                   by, n), call. = FALSE)
    }
    facet <- n > 4
    pal <- if (facet) rep(maf_colours(1), n) else colours %||% maf_colours(n)
    if (length(pal) < n) stop("`colours` needs ", n, " colours.", call. = FALSE)
    pal <- stats::setNames(pal[seq_len(n)], levels(cls))
    if (anyNA(cls)) {
      cls <- factor(cls, levels = c(levels(cls), "Unknown"))
      cls[is.na(cls)] <- "Unknown"
      pal <- c(pal, Unknown = if (facet) maf_colours(1) else .ink$unknown)
    }
  }

  # Expected -log10(P) of the k-th smallest of n uniform P-values, within each
  # class: sort once by class, then by P.
  g <- as.integer(cls)
  n_k <- tabulate(g, nlevels(cls))
  o <- order(g, -logp, method = "radix")
  lp <- logp[o]
  g <- g[o]
  start <- c(0, cumsum(n_k))
  expected <- -log10((seq_along(lp) - start[g] - 0.5) / n_k[g])

  # Genomic inflation, from the median P-value: overall, and per class from the
  # middle of each class's sorted block.
  lam <- .lambda_gc(stats::median(logp))
  lam_k <- vapply(seq_along(n_k), function(j) {
    n <- n_k[j]
    if (!n) return(NA_real_)
    s <- start[j]
    .lambda_gc(if (n %% 2) lp[s + (n + 1) / 2] else (lp[s + n / 2] + lp[s + n / 2 + 1]) / 2)
  }, 0)

  max_logp <- max(lp)
  comp <- .y_compressor(max_logp, sig, y_compress)
  keep <- .thin(expected, lp, below = sig, group = g)
  pd <- data.frame(x = expected[keep], y = comp$f(lp[keep]),
                   cls = factor(levels(cls)[g[keep]], levels(cls)), g = g[keep])
  # Rarest drawn last (on top); unknown frequency first (underneath).
  pd <- pd[order(ifelse(pd$cls %in% "Unknown", 0L, pd$g)), ]

  present <- levels(cls)[n_k > 0]
  bands <- do.call(rbind, lapply(which(n_k > 0), function(j) {
    k <- unique(round(10^seq(0, log10(n_k[j]), length.out = 500)))
    data.frame(cls = levels(cls)[j], x = -log10((k - 0.5) / n_k[j]),
               lo = comp$f(-log10(stats::qbeta(0.975, k, n_k[j] - k + 1))),
               hi = comp$f(-log10(stats::qbeta(0.025, k, n_k[j] - k + 1))))
  }))
  xmax <- max(expected)
  diagonal <- data.frame(x = seq(0, xmax, length.out = 200))
  diagonal$y <- comp$f(diagonal$x)
  ymax <- max(comp$f(max_logp), comp$f(sig), xmax) * 1.06
  yb <- .y_breaks(comp, max_logp)
  pts <- point_size * base_size / 10

  fmt <- function(v) formatC(v, format = "f", digits = 3)
  # Per-class lambda only where it is precise enough to read (SE ~ 2.3/sqrt(n)):
  # in the legend, or in each panel's title.
  show_lam <- lambda & n_k >= 1000
  leg_lab <- .maf_legend_labels(present)
  if (any(show_lam)) {
    leg_lab <- as.expression(lapply(seq_along(present), function(i) {
      j <- match(present[i], levels(cls))
      if (!show_lam[j]) leg_lab[[i]] else bquote(.(leg_lab[[i]]) ~~~ lambda == .(fmt(lam_k[j])))
    }))
  }
  if (facet) {
    strip <- vapply(present, function(l) {
      j <- match(l, levels(cls))
      paste0(.pm_string(l), if (show_lam[j]) paste0(" ~~~ lambda == '", fmt(lam_k[j]), "'"))
    }, "")
    pd$strip <- factor(strip[as.character(pd$cls)], levels = unname(strip))
    bands$strip <- factor(strip[as.character(bands$cls)], levels = unname(strip))
  }

  p <- ggplot2::ggplot()
  if (ci) {
    p <- p + ggplot2::geom_ribbon(
      data = bands, ggplot2::aes(x = .data$x, ymin = .data$lo, ymax = .data$hi, group = .data$cls),
      fill = .ink$noise, alpha = 0.2, colour = NA)
  }
  if (!is.null(comp$at)) p <- p + .axis_break(comp$at, ymax)
  p <- p +
    ggplot2::geom_line(data = diagonal, ggplot2::aes(x = .data$x, y = .data$y),
                       colour = .ink$secondary, linewidth = 0.35) +
    ggplot2::geom_hline(yintercept = sig, colour = .ink$secondary, linewidth = 0.3,
                        linetype = "22") +
    ggplot2::geom_point(data = pd, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$cls),
                        size = 0.9 * pts, stroke = 0) +
    ggplot2::scale_colour_manual(
      values = pal, limits = present, labels = leg_lab,
      name = if (by == "maf") "Minor allele frequency" else by,
      guide = if (by != "none" && !facet) {
        ggplot2::guide_legend(nrow = if (length(present) > 2) 2 else 1, byrow = TRUE,
                              override.aes = list(size = 2.6 * pts))
      } else "none"
    )
  if (facet) {
    p <- p + ggplot2::facet_wrap(~strip, labeller = ggplot2::label_parsed)
    if (lambda) p <- p + ggplot2::labs(subtitle = bquote(lambda[GC] == .(fmt(lam)) ~ "for all variants"))
  } else if (lambda) {
    p <- p + ggplot2::annotate("text", x = 0, y = ymax, hjust = 0, vjust = 1, parse = TRUE,
                               label = sprintf("lambda[GC] == '%s'", fmt(lam)),
                               size = base_size * .text_size$legend / ggplot2::.pt, colour = .ink$primary)
  }
  p +
    ggplot2::scale_x_continuous(limits = c(0, xmax * 1.02), expand = ggplot2::expansion(0)) +
    ggplot2::scale_y_continuous(
      breaks = yb$at, labels = yb$labels, limits = c(0, ymax),
      expand = ggplot2::expansion(0),
      sec.axis = ggplot2::dup_axis(name = NULL, breaks = sig,
                                   labels = parse(text = .p_plotmath(threshold)))
    ) +
    ggplot2::labs(x = expression(Expected ~ -log[10](italic(P))),
                  y = expression(Observed ~ -log[10](italic(P)))) +
    theme_brownstone(base_size = base_size) + .legend_on_top(base_size) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(colour = .ink$grid, linewidth = 0.3),
      axis.text.y.right = ggplot2::element_text(colour = .ink$muted, size = base_size * .text_size$note,
                                                hjust = 0)
    )
}

# Genomic inflation from the median -log10(P): the 1-df chi-square at the
# median P-value over its expected median, on the log scale so that it never
# underflows.
.lambda_gc <- function(median_logp) {
  stats::qchisq(-median_logp * log(10), df = 1, lower.tail = FALSE, log.p = TRUE) /
    stats::qchisq(0.5, df = 1)
}
