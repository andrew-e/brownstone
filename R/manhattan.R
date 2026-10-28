#' Manhattan plot with effect sizes, allele frequency and genes
#'
#' @description
#' Draws a Manhattan plot that answers three questions a classic one can't.
#'
#' * **How big is each effect, and how precise?** Beneath the chromosome axis an
#'   aligned track shows the size of each lead variant's effect (|beta|, or the
#'   odds ratio of the risk allele) with its 95% confidence interval, so every
#'   tower of P-values stands directly above its own effect estimate.
#' * **Is the signal from a common or a rare variant?** Significant variants are
#'   coloured by minor allele frequency, by class (blue, aqua, orange) or on a
#'   continuous scale. The effect track also shows the strongest variant of each frequency
#'   class within a locus, whenever it is too significant to be an LD shadow of
#'   the locus lead, so allelic series (a common and a rare signal at the same
#'   gene) become visible.
#' * **Which genes?** Each locus is labelled with its nearest gene, or with your
#'   own labels.
#'
#' @param data GWAS summary statistics, one row per variant. Columns are
#'   detected automatically; see [standardise_sumstats()].
#' @param ... Column names and options passed to [standardise_sumstats()], for
#'   example `chr = "CHROM", pos = "GENPOS", eaf = "A1FREQ"`, or
#'   `effect_allele = "as_is"`.
#' @param build Genome build of the positions, used for gene labels.
#' @param threshold Genome-wide significance threshold.
#' @param suggestive Optional suggestive threshold (e.g. `1e-5`). Adds a faint
#'   line, draws peaks that reach it but not `threshold` as hollow rings in
#'   their frequency colour, labels them in grey (`label_suggestive`) and puts
#'   them on the effect track as hollow points (`track_suggestive`). Variants
#'   between the lines that belong to a genome-wide significant locus stay
#'   grey: they are the base of that locus's tower, not a peak of their own.
#' @param label_suggestive Label suggestive loci (after all genome-wide
#'   significant ones, within `max_labels`).
#' @param track_suggestive Show suggestive loci on the effect track. Their
#'   effects are the most inflated by winner's curse, and a rare variant that
#'   only reaches the suggestive line sits high on the track simply because
#'   its standard error is large.
#' @param labels How to label loci: `"nearest_gene"` (the default), `"snp"`
#'   (variant IDs), the name of a column in `data` (its value at the lead
#'   variant), a data frame of your own labels, or `NULL` for none. A label data
#'   frame needs a `label` column and either `snp` or `chr` + `pos`; an optional
#'   `fontface` column (`"plain"`, `"bold"`, `"italic"`, `"bold.italic"`) styles
#'   individual labels, e.g. to set novel loci in bold.
#' @param max_labels Label at most this many loci (the most significant).
#' @param window Half-width, in base pairs, used to group significant variants
#'   into loci; see [find_loci()].
#' @param effects What the effect track shows (it needs effect and standard
#'   error columns). `"magnitude"` (the default): |beta| with its 95% CI, i.e.
#'   the effect per copy of the trait-increasing allele, or for odds ratios
#'   the OR of the risk allele. The sign of a single variant's effect only
#'   reflects which allele the GWAS reported, so at a glance the size is what
#'   matters. A lead variant whose effect disagrees with its own case/control
#'   allele frequencies (by an order of magnitude, or by six standard errors;
#'   binary traits whose summary statistics report them) is marked with a
#'   black ring; see [find_artifacts()].
#'   `"signed"`: beta with its sign, oriented to the minor allele when an allele
#'   frequency is available (see `effect_allele` in [standardise_sumstats()]),
#'   for lead variants only (`track_variants` is ignored: LD proxies are
#'   reported against arbitrary alleles, so their signs mean nothing next to
#'   the lead's).
#'   `"none"` (or `FALSE`): no track.
#' @param effect_scale How the effect track's axis reacts to the data.
#'   `"auto"` (the default): linear, so effects are drawn in proportion, unless
#'   a few estimates dwarf the rest. Then, as on the P-value axis, the axis
#'   breaks (at 3 times the median lead effect, rounded) and the outliers go on
#'   a log scale in the top 25-45% of the track. The break is triggered by more
#'   than twice that value among the point estimates; wide confidence
#'   intervals alone don't count. `"sqrt"`: a smooth square-root scale.
#'   `"linear"`: never break.
#' @param allelic_series If `TRUE` (the default), the effect track adds the
#'   strongest variant of each other frequency class in a locus, when its
#'   z-statistic exceeds the most that LD with the locus lead could produce
#'   (allele frequencies cap the possible LD). If `FALSE`, the track shows one
#'   variant per locus.
#' @param maf_breaks Boundaries of the minor allele frequency classes. The
#'   default, `c(0.001, 0.01, 0.05)`, gives ultra-rare (< 0.1%), rare
#'   (0.1-1%), low-frequency (1-5%) and common (>= 5%) classes.
#' @param maf_scale `"binned"` (the default) colours significant variants by
#'   frequency class (see `maf_breaks`), one distinct hue per class, so that
#'   rare variants stand out. `"continuous"` colours them by log10(MAF) on a
#'   colour ramp, with no cut-offs: from the rarest hit's order of magnitude
#'   (at least down to 0.1%) up to 10%, above which all common variants share
#'   one colour (differences between 10% and 40% don't matter here, and
#'   spending colour on them only mottles the towers). Either way the colour
#'   follows the minor allele frequency, which, unlike the effect allele
#'   frequency, doesn't depend on which allele the GWAS reported. Both default
#'   palettes are colour-blind safe.
#' @param colours For `maf_scale = "binned"`, one colour per frequency class,
#'   most common first (default [maf_colours()]); for `"continuous"`, a colour
#'   ramp from common to rare (default [maf_ramp()]).
#' @param y_compress Where to break the y axis, in -log10(P). Below the break
#'   the axis is linear; above it -log10(P) is on a log scale given about 15%
#'   of the panel height per decade (12-35% in total, and never more than a
#'   linear axis would give it), so the long tails keep their order without
#'   taking space from the threshold region. `"auto"` (the default) breaks at
#'   2.5 times the significance line, rounded (20 for 5e-8), and only once the
#'   data climb to about twice that line; a number sets the break; `FALSE`
#'   turns it off.
#' @param track_height Fraction of the figure height given to the effect track.
#' @param track_variants `"all"` (the default) shows every variant drawn in
#'   colour in the main panel: lead variants as large points with their 95%
#'   CI, every other hit as a small point without a CI (dots for genome-wide
#'   significant variants, rings for suggestive peaks, grey for the bases of
#'   towers). Most of these are LD proxies, whose expected |beta| is the
#'   lead's shrunk by their correlation, so they sit just below the lead. LD
#'   with the lead can never lift a proxy's expected |beta| above the lead's,
#'   so a point clearly above its lead (allowing for the wide CIs of rare
#'   variants) points to a second signal. Caution: rare variants that only
#'   reach the suggestive line sit high simply because their standard errors
#'   are large. `"leads"` shows only the lead variants (one number per locus,
#'   plus any allelic-series extras), and is always used for signed effects.
#' @param track_direction `"up"` (the default) draws the effect track rising
#'   from its own baseline; `"down"` hangs it from the Manhattan's axis, like a
#'   reflection of the skyline.
#' @param label_angle `0` for horizontal labels, `90` for vertical labels
#'   (which pack better when there are many loci).
#' @param base_size Base font size, in points.
#' @param point_size Multiplier for point sizes.
#'
#' @return A patchwork object, or a ggplot when there is no effect track. Add a
#'   title with `+ patchwork::plot_annotation(title = "...")`.
#' @export
#' @examples
#' gwas <- simulate_gwas(n_variants = 2e5)
#' manhattan(gwas)
manhattan <- function(data, ..., build = c("GRCh38", "GRCh37"),
                          threshold = 5e-8, suggestive = NULL,
                          labels = "nearest_gene", label_suggestive = TRUE,
                          track_suggestive = TRUE,
                          max_labels = 40, window = 5e5,
                          effects = c("magnitude", "signed", "none"),
                          effect_scale = c("auto", "sqrt", "linear"),
                          allelic_series = TRUE,
                          maf_scale = c("binned", "continuous"),
                          maf_breaks = c(0.001, 0.01, 0.05), colours = NULL,
                          y_compress = "auto", track_height = 0.22,
                          track_variants = c("all", "leads"),
                          track_direction = c("up", "down"),
                          label_angle = 0, base_size = 10, point_size = 1) {
  build <- match.arg(build)
  effect_scale <- match.arg(effect_scale)
  track_direction <- match.arg(track_direction)
  track_variants <- match.arg(track_variants)
  maf_scale <- match.arg(maf_scale)
  if (isTRUE(effects)) effects <- "magnitude"
  if (is.null(effects) || isFALSE(effects)) effects <- "none"
  effects <- match.arg(effects)
  if (effects == "signed" && track_variants == "all") {
    # Proxies are reported against arbitrary alleles, so their signs are
    # meaningless next to the lead's: signed effects show lead variants only.
    track_variants <- "leads"
  }
  label_col <- if (is.character(labels) && length(labels) == 1L &&
                    !labels %in% c("nearest_gene", "snp")) labels
  d <- .standardise(data, ..., keep = label_col, ids = FALSE, ccafs = TRUE)
  sig <- -log10(threshold)
  if (!is.null(suggestive) && suggestive <= threshold) {
    stop("`suggestive` must be a larger P-value than `threshold`.", call. = FALSE)
  }
  sugg <- if (!is.null(suggestive)) -log10(suggestive)
  floor_logp <- sugg %||% sig   # variants at or above this are "hits"

  layout <- .genome_layout(d)
  d$x <- .genome_x(d$chr, d$pos, layout)
  max_logp <- max(d$logp)
  # A table of labels can name any variant, so it is matched before thinning.
  if (is.data.frame(labels)) labels <- .label_table(labels, d)
  # Below the lowest line only one point per cell of a fine grid is drawn (see
  # .thin()). Every hit is kept, so all that follows works on far fewer rows.
  d <- d[.thin(d$x, d$logp, below = floor_logp), ]
  # Lead variants whose effects disagree with their own case/control allele
  # frequencies (see find_artifacts()), marked with a ring on the effect track.
  flag_rows <- .artifact_flags(d, threshold, window, build)
  has_maf <- isTRUE(attr(d, "has_maf"))
  cset <- .colour_setup(d, has_maf, maf_scale, maf_breaks, colours, floor_logp)
  cs <- cset$cs
  d[c("class", "col", "stack")] <- cset[c("class", "col", "stack")]

  loc <- .locus_setup(d, floor_logp, sig, window, labels, label_col, label_suggestive,
                      build, max_labels)
  d <- loc$d
  cl <- loc$cl
  lab <- loc$lab
  comp <- .y_compressor(max_logp, sig, y_compress)

  on_track <- loc$strong_locus | (track_suggestive & !is.null(sugg))
  has_effects <- effects != "none" && isTRUE(attr(d, "has_effect")) && any(on_track)
  if (effects != "none" && !isTRUE(attr(d, "has_effect"))) {
    message("No effect size and standard error columns found; drawing without the effect track.")
  }

  # Over a rising track, the floor stands a little above the ground line
  # (.ground_rule()), so that the line reads as its own; a hanging track hangs
  # straight from it.
  main <- .main_panel(d, max_logp, layout, comp, sig, sugg, lab, cs, has_maf,
                      label_angle, base_size, point_size,
                      x_axis = !has_effects,
                      floor_gap = if (has_effects && track_direction == "up") 0.03 else 0)
  if (!has_effects) return(main)

  tr <- .track_setup(d, cl, sig, floor_logp, window, allelic_series && has_maf,
                     track_suggestive, track_variants, flag = flag_rows)
  eff <- tr$eff
  cloud <- tr$cloud
  track <- .effect_panel(eff, cloud, layout, cs, mode = effects, scale = effect_scale,
                         log_or = identical(attr(d, "effect_scale"), "log_or"),
                         allele = attr(d, "effect_allele"),
                         down = track_direction == "down",
                         base_size = base_size, point_size = point_size)
  # A bold ground line under the towers, and a little air before a track that
  # rises from its own baseline (a hanging track hangs from that line).
  main <- main + .ground_rule() + ggplot2::theme(plot.margin = ggplot2::margin(6, 10, 0, 6))
  gap <- if (track_direction == "down") 2 else 8
  track <- track + ggplot2::theme(plot.margin = ggplot2::margin(gap, 10, 4, 6))
  # Titles added with plot_annotation() take the theme's (large, centred) style.
  patchwork::wrap_plots(main, track, ncol = 1,
                        heights = c(1 - track_height, track_height)) +
    patchwork::plot_annotation(theme = theme_maffect(base_size = base_size))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Frequency classes, palette and colour scales for the (thinned) variants in
# `d`: the class of each row, the value its colour is mapped from (`col`), its
# drawing order (`stack`) and the scales (`cs`, see .colour_spec()).
# miami_plot() calls this once on both GWAS together, so that the two halves
# share one palette.
.colour_setup <- function(d, has_maf, maf_scale, maf_breaks, colours, floor_logp) {
  cls <- if (has_maf) .maf_classes(d$maf, maf_breaks) else factor(rep("all", nrow(d)))
  n_cls <- nlevels(cls)
  binned_cols <- if (maf_scale == "binned") colours
  pal <- if (has_maf) binned_cols %||% maf_colours(n_cls) else (binned_cols %||% maf_colours(1))[1]
  if (length(pal) < n_cls) {
    stop("`colours` needs ", n_cls, " colours.", call. = FALSE)
  }
  pal <- stats::setNames(pal[seq_len(n_cls)], levels(cls))
  # Significant variants with no frequency get their own neutral class rather
  # than silently looking like noise.
  if (has_maf && anyNA(cls[d$logp >= floor_logp])) {
    cls <- factor(cls, levels = c(levels(cls), "Unknown"))
    cls[is.na(cls)] <- "Unknown"
    pal <- c(pal, Unknown = .ink$unknown)
  }
  d$class <- cls
  cs <- .colour_spec(d, has_maf, maf_scale, pal,
                     if (maf_scale == "continuous") colours, floor_logp)
  list(class = cls, col = cs$values, stack = cs$stack, cs = cs)
}

# Loci (found down to the lowest line) and their labels. Adds
# `in_strong_locus` to `d`.
.locus_setup <- function(d, floor_logp, sig, window, labels, label_col, label_suggestive,
                         build, max_labels) {
  # Genome-wide significant leads come first because clumping works in order
  # of significance.
  cl <- .clump(d, floor_logp, window)
  # Variants between the two lines that belong to a genome-wide significant
  # locus are just the base of its tower, not a suggestive peak.
  strong_locus <- d$logp[cl$lead] >= sig
  d$in_strong_locus <- FALSE
  d$in_strong_locus[cl$member[strong_locus[cl$locus_of]]] <- TRUE
  label_leads <- if (label_suggestive) cl$lead else cl$lead[strong_locus]
  lab <- .locus_labels(labels, label_col, d, label_leads, sig, build, max_labels)
  list(d = d, cl = cl, strong_locus = strong_locus, lab = lab)
}

# Lead variants (`eff`) and, for `track_variants = "all"`, every other hit
# (`cloud`) for the effect track. `flag`: the rows of `d` whose lead variants
# disagree with their own case/control allele frequencies (see
# find_artifacts()), marked with a ring on the track.
.track_setup <- function(d, cl, sig, floor_logp, window, allelic_series, track_suggestive,
                         track_variants, flag = integer(0)) {
  eff <- .track_variants(d, cl, sig, window, allelic_series, suggestive = track_suggestive)
  eff$flag <- eff$row %in% flag
  cloud <- if (track_variants == "all") {
    .track_cloud(d, eff$row, if (track_suggestive) floor_logp else sig, sig)
  }
  list(eff = eff, cloud = cloud)
}

# Minor allele frequency classes, as a factor ordered from common to rare.
# Levels are ASCII ("< 1%", "1-5%", ">= 5%"); .maf_legend_labels() turns them
# into typeset legend labels.
.maf_classes <- function(maf, breaks) {
  breaks <- sort(unique(breaks))
  num <- function(b) format(b * 100, drop0trailing = TRUE, trim = TRUE)
  n <- length(breaks)
  labs <- c(paste0("< ", num(breaks[1]), "%"),
            if (n > 1) paste0(num(breaks[-n]), "-", num(breaks[-1]), "%"),
            paste0(">= ", num(breaks[n]), "%"))
  f <- cut(maf, c(0, breaks, 0.5 + 1e-9), labels = labs, right = FALSE,
           include.lowest = TRUE)
  # Reverse the levels by reversing the codes: factor(f, levels = rev(labs))
  # gives the same, but goes through text, which is slow on millions of rows.
  structure(length(labs) + 1L - as.integer(f), levels = rev(labs), class = "factor")
}

# Legend labels via plotmath, so that ">=" is set as a proper symbol on every
# graphics device (raw Unicode text fails on the default pdf() device).
.maf_legend_labels <- function(levels) {
  as.expression(lapply(levels, function(l) {
    if (startsWith(l, ">= ")) bquote("" >= .(substring(l, 4))) else l
  }))
}

# How points are coloured, shared by both panels: the column to map (`values`),
# the drawing order (`stack`: rarest drawn last, on top) and the fill and colour
# scales. Binned: one colour per frequency class. Continuous: log10(MAF) on a
# colour ramp, with conventional cut-offs as ticks.
.colour_spec <- function(d, has_maf, maf_scale, pal, ramp, floor_logp) {
  if (!has_maf || maf_scale == "binned") {
    fill <- function(legend = TRUE) {
      ggplot2::scale_fill_manual(
        values = pal, labels = .maf_legend_labels(names(pal)), drop = FALSE,
        name = "Minor allele frequency",
        guide = if (legend) {
          ggplot2::guide_legend(order = 1, override.aes = list(shape = 21, colour = "white",
                                                               size = 3))
        } else "none"
      )
    }
    colour <- function(legend = FALSE) {
      ggplot2::scale_colour_manual(
        values = pal, labels = .maf_legend_labels(names(pal)), drop = FALSE,
        name = "Minor allele frequency",
        guide = if (legend) {
          ggplot2::guide_legend(order = 1, override.aes = list(shape = 21, fill = "white",
                                                               stroke = 0.8, size = 2.6))
        } else "none"
      )
    }
    # Rarest drawn last (on top); unknown frequency first (underneath).
    stack <- ifelse(d$class %in% "Unknown", 0L, as.integer(d$class))
    return(list(values = d$class, stack = stack, fill = fill, colour = colour))
  }
  ramp <- ramp %||% maf_ramp()
  # The rare end of the scale follows the data: down to the rarest hit's
  # order of magnitude, and at least down to 0.1%.
  hit_maf <- d$maf[d$logp >= floor_logp & is.finite(d$maf) & d$maf > 0]
  lo <- min(1e-3, if (length(hit_maf)) 10^floor(log10(min(hit_maf))) else 1e-3)
  # Above 10% the colour stops changing: differences between common variants
  # (10% vs 40%) don't matter here, and spending colour on them only mottles
  # the towers. The gradient goes where frequency matters, rare to common.
  hi <- 0.1
  lim <- log10(c(lo, hi))
  brk <- c(1e-5, 1e-4, 1e-3, 1e-2, hi)
  brk <- brk[brk >= lo]
  # Labels go through plotmath so ">=" is a real symbol on every device.
  pct <- paste0(format(brk * 100, drop0trailing = TRUE, trim = TRUE,
                       scientific = FALSE), "%")
  lab <- .as_plotmath(c(sprintf("'%s'", pct[-length(pct)]),
                        sprintf("'' >= '%s'", pct[length(pct)])))
  gradient <- function(f, legend) {
    f(colours = rev(ramp), limits = lim, breaks = log10(brk), labels = lab,
      oob = scales::oob_squish, na.value = .ink$unknown,
      name = "Minor allele frequency",
      guide = if (legend) {
        ggplot2::guide_colourbar(order = 1, theme = ggplot2::theme(
          legend.key.width = ggplot2::unit(8, "lines"),
          legend.key.height = ggplot2::unit(0.45, "lines"),
          legend.title.position = "left",
          legend.title = ggplot2::element_text(vjust = 1, margin = ggplot2::margin(r = 6))))
      } else "none")
  }
  list(values = log10(pmin(pmax(d$maf, lo), hi)),
       stack = ifelse(is.na(d$maf), -Inf, -d$maf),
       fill = function(legend = TRUE) gradient(ggplot2::scale_fill_gradientn, legend),
       colour = function(legend = FALSE) gradient(ggplot2::scale_colour_gradientn, legend))
}

# Distance clumping on a standardised table. Returns the row index of each
# lead variant and, for every significant variant, the locus it belongs to.
.clump <- function(d, min_logp, window) {
  sig <- which(d$logp >= min_logp)
  o <- sig[order(-d$logp[sig])]
  chr_o <- d$chr[o]
  pos_o <- d$pos[o]
  locus_of <- rep(NA_integer_, length(o))
  leads <- integer(0)
  for (i in seq_along(o)) {
    if (!is.na(locus_of[i])) next
    leads <- c(leads, i)
    claim <- is.na(locus_of) & chr_o == chr_o[i] & abs(pos_o - pos_o[i]) <= window
    locus_of[claim] <- length(leads)
  }
  list(lead = o[leads], member = o, locus_of = locus_of)
}

# |z| from -log10(P), without underflow.
.z_from_logp <- function(logp) {
  abs(stats::qnorm(-logp * log(10) - log(2), log.p = TRUE))
}

# The largest possible correlation between two biallelic variants with minor
# allele frequencies f1 and f2.
.r_max <- function(f1, f2) {
  lo <- pmin(f1, f2)
  hi <- pmax(f1, f2)
  sqrt(lo * (1 - hi) / (hi * (1 - lo)))
}

# Variants for the effect track: the lead of every genome-wide significant
# locus (and, with `suggestive = TRUE`, of every suggestive locus), plus (for
# allelic series) the lead of each other frequency class in a genome-wide
# significant locus when it is too significant to be explained by LD with the
# locus lead.
.track_variants <- function(d, cl, sig, window, allelic_series, suggestive = FALSE) {
  idx <- if (suggestive) cl$lead else cl$lead[d$logp[cl$lead] >= sig]
  if (allelic_series && length(cl$member)) {
    locus <- cl$locus_of
    for (lev in setdiff(levels(d$class), "Unknown")) {
      in_class <- cl$member[d$class[cl$member] %in% lev]
      if (!length(in_class)) next
      sub <- .clump(d[in_class, c("chr", "pos", "logp")], sig, window)
      cand <- in_class[sub$lead]
      lead_of_cand <- cl$lead[locus[match(cand, cl$member)]]
      cand_ok <- cand != lead_of_cand &
        .z_from_logp(d$logp[cand]) >
          .r_max(d$maf[cand], d$maf[lead_of_cand]) * .z_from_logp(d$logp[lead_of_cand]) + 2
      idx <- c(idx, cand[!is.na(cand_ok) & cand_ok])
    }
  }
  idx <- unique(idx)
  out <- d[idx, c("chr", "pos", "x", "snp", "logp", "beta", "se", "maf", "class",
                  "col", "stack")]
  out$row <- idx
  out$lead <- idx %in% cl$lead
  out$strong <- out$logp >= sig
  out[is.finite(out$beta) & is.finite(out$se), ]
}

# Every other variant drawn in colour (or as a ring, or as the grey base of a
# tower) in the main panel, for `track_variants = "all"`: the same points get
# the same look in both panels. Thinned to one point per cell of a fine grid.
.track_cloud <- function(d, lead_rows, floor_logp, sig) {
  i <- which(d$logp >= floor_logp & is.finite(d$beta) & is.finite(d$se))
  i <- setdiff(i, lead_rows)
  out <- d[i, c("x", "beta", "se", "class", "col", "stack", "logp", "in_strong_locus")]
  out$tier <- factor(ifelse(out$logp >= sig, "gws",
                            ifelse(out$in_strong_locus, "base", "suggestive")),
                     levels = c("base", "suggestive", "gws"))
  if (!nrow(out)) return(out)
  keep <- .thin(out$x, sqrt(abs(out$beta)), below = Inf, nx = 3000, ny_per_unit = 300,
                group = (as.integer(out$class) - 1L) * 3L + as.integer(out$tier))
  out[keep, ]
}

# Label table for the main panel. Labels of loci that are only suggestive are
# set in a quieter ink than genome-wide significant ones.
.locus_labels <- function(labels, label_col, d, leads, sig, build, max_labels) {
  out <- .locus_label_text(labels, label_col, d, leads, build, max_labels)
  out$colour <- ifelse(out$logp >= sig, .ink$primary, .ink$muted)
  out
}

.locus_label_text <- function(labels, label_col, d, leads, build, max_labels) {
  empty <- data.frame(x = numeric(0), logp = numeric(0), label = character(0),
                      fontface = character(0), stringsAsFactors = FALSE)
  if (is.null(labels) || isFALSE(labels)) return(empty)
  if (is.data.frame(labels)) return(labels)   # already matched by .label_table()
  lead <- utils::head(leads, max_labels)
  if (!length(lead)) return(empty)
  if (identical(labels, "nearest_gene")) {
    txt <- nearest_gene(d$chr[lead], d$pos[lead], build = build)$gene
    txt[is.na(txt)] <- .variant_id(d, lead)[is.na(txt)]
    ff <- "italic"
  } else if (identical(labels, "snp")) {
    txt <- .variant_id(d, lead)
    ff <- "plain"
  } else {
    txt <- as.character(d[[label_col]][lead])
    ff <- "plain"
  }
  data.frame(x = d$x[lead], logp = d$logp[lead], label = txt, fontface = ff,
             stringsAsFactors = FALSE)
}

# A table of labels, matched to the data by variant ID or by position. Any
# variant can be labelled, so this runs on every row, before thinning.
.label_table <- function(labels, d, warn = TRUE) {
  if (!"label" %in% names(labels)) {
    stop("A label data frame needs a `label` column.", call. = FALSE)
  }
  if (!"snp" %in% names(labels) && !all(c("chr", "pos") %in% names(labels))) {
    stop("A label data frame needs a `snp` column, or `chr` and `pos` columns.",
         call. = FALSE)
  }
  i <- .match_variants(labels, d)
  if (warn && anyNA(i)) {
    warning(sum(is.na(i)), " labelled variant(s) not found in the data.", call. = FALSE)
  }
  ff <- if ("fontface" %in% names(labels)) as.character(labels$fontface) else "plain"
  out <- data.frame(x = d$x[i], logp = d$logp[i], label = as.character(labels$label),
                    fontface = rep_len(ff, length(i)), stringsAsFactors = FALSE)
  out[!is.na(i), ]
}

# Chromosome numbers at the note size, the smallest text of the package: the
# short chromosomes 19-22 are close together, and at the axis size a wider
# typeface (Avenir Next) loses 21 to the overlap check at 8 in.
.chr_axis_text <- function(base_size) {
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = base_size * .text_size$note))
}

# `top_title` adds a title above the panel, below any legend (miami_plot()
# names its upper half this way).
.x_scale <- function(layout, labels = TRUE, top_title = NULL) {
  ggplot2::scale_x_continuous(
    breaks = layout$mid,
    labels = if (labels) layout$chr else NULL,
    limits = c(min(layout$start), max(layout$end)),
    expand = ggplot2::expansion(0),
    guide = ggplot2::guide_axis(check.overlap = TRUE),
    sec.axis = if (is.null(top_title)) ggplot2::waiver() else {
      ggplot2::dup_axis(name = top_title, breaks = NULL, labels = NULL)
    }
  )
}

.band_layer <- function(layout) {
  b <- layout[layout$index %% 2 == 0, ]
  ggplot2::annotate("rect", xmin = b$start, xmax = b$end, ymin = -Inf, ymax = Inf,
                    fill = .ink$band, colour = NA)
}

# Break the panel where the -log10(P) axis switches to its compressed scale
# (see .y_compressor()), with the broken-axis mark (see .break_mark()).
.axis_break <- function(at, ymax) {
  gap <- 0.008 * ymax
  .break_mark(at - gap, at + gap, "y")
}

# "5e-08" -> plotmath for 5 x 10^-8
.p_plotmath <- function(p) {
  e <- floor(log10(p))
  m <- signif(p / 10^e, 2)
  if (m == 1) sprintf("10^%d", e) else sprintf("%s %%*%% 10^%d", format(m), e)
}

# `d` is already thinned (see manhattan()); `max_logp` is the top before thinning.
# With `down = TRUE` the panel hangs from its top edge (the lower half of a
# Miami plot): -log10(P) is drawn negated and labelled as positive, which
# keeps ggrepel's nudges and limits simple.
.main_panel <- function(d, max_logp, layout, comp, sig, sugg, lab, cs, has_maf,
                        label_angle, base_size, point_size, x_axis, floor_gap = 0,
                        down = FALSE, legend = TRUE, top_title = NULL) {
  s <- if (down) -1 else 1
  floor_logp <- sugg %||% sig
  pd <- d[, c("x", "logp", "class", "col", "stack", "in_strong_locus")]
  pd$y <- s * comp$f(pd$logp)
  tower_base <- pd$logp >= floor_logp & pd$logp < sig & pd$in_strong_locus
  noise <- pd[pd$logp < floor_logp | tower_base, ]
  hits <- pd[pd$logp >= floor_logp & !tower_base, ]
  hits$tier <- factor(ifelse(hits$logp >= sig, "gws", "suggestive"),
                      levels = c("gws", "suggestive"))
  hits <- hits[order(hits$stack), ]   # rarest last, i.e. on top
  strong <- hits[hits$tier == "gws", ]
  weak <- hits[hits$tier == "suggestive", ]

  ytop <- comp$f(max_logp)
  ymax <- max(ytop * 1.1, sig * 1.6)
  yb <- .y_breaks(comp, max_logp)
  pts <- point_size * base_size / 10

  p <- ggplot2::ggplot() + .band_layer(layout)
  if (!is.null(comp$at)) p <- p + .axis_break(s * comp$at, ymax)
  p <- p +
    ggplot2::geom_point(data = noise, ggplot2::aes(x = .data$x, y = .data$y),
                        colour = .ink$noise, size = 0.45 * pts, stroke = 0)
  if (!is.null(sugg)) {
    p <- p + ggplot2::geom_hline(yintercept = s * sugg, colour = .ink$muted,
                                 linewidth = 0.25, linetype = "13")
  }
  # Suggestive hits are hollow rings in their frequency colour; genome-wide
  # significant hits are solid. Both map `shape` to the tier only to get a
  # legend entry for it.
  p <- p +
    ggplot2::geom_hline(yintercept = s * sig, colour = .ink$secondary, linewidth = 0.3,
                        linetype = "22")
  # With no hits at all there is nothing to colour (and empty manual scales
  # only produce warnings). Legend keys are drawn for every class, even one
  # with no hits in this panel (the other half of a Miami plot may have them).
  if (nrow(hits)) p <- p +
    ggplot2::geom_point(data = weak, ggplot2::aes(x = .data$x, y = .data$y,
                                                  colour = .data$col, shape = .data$tier),
                        fill = "white", stroke = 0.6 * pts, size = 1.45 * pts,
                        show.legend = TRUE) +
    ggplot2::geom_point(data = strong, ggplot2::aes(x = .data$x, y = .data$y,
                                                    fill = .data$col, shape = .data$tier),
                        colour = "white", stroke = 0.25 * pts, size = 1.7 * pts,
                        show.legend = TRUE) +
    # The frequency legend comes from the solid points, or from the rings when
    # nothing is genome-wide significant.
    (if (nrow(strong)) cs$fill(legend = legend && has_maf)) +
    cs$colour(legend = legend && has_maf && !nrow(strong)) +
    ggplot2::scale_shape_manual(
      values = c(gws = 21, suggestive = 21), limits = c("gws", "suggestive"),
      labels = parse(text = paste("italic(P) <", vapply(10^-c(sig, sugg %||% sig), .p_plotmath, ""))),
      name = NULL,
      guide = if (legend && !is.null(sugg)) {
        ggplot2::guide_legend(order = 2, override.aes = list(
          fill = c(.ink$secondary, "white"), colour = c("white", .ink$secondary),
          stroke = c(0.25, 0.6), size = c(3, 2.6)))
      } else "none"
    )

  if (nrow(lab)) {
    lab$y <- s * comp$f(lab$logp)
    # Labels repel from every significant point, not just the labelled ones;
    # a coarse grid of hits is enough for that and keeps ggrepel fast.
    near <- hits[.thin(hits$x, hits$y, below = Inf, nx = 600, ny_per_unit = 1), ]
    rep <- rbind(lab[, c("x", "y", "label", "fontface", "colour")],
                 data.frame(x = near$x, y = near$y, label = "", fontface = "plain",
                            colour = .ink$primary))
    common <- list(
      data = rep, mapping = ggplot2::aes(x = .data$x, y = .data$y, label = .data$label,
                                         fontface = .data$fontface),
      colour = rep$colour, size = base_size * .text_size$label / ggplot2::.pt,
      family = .family(), min.segment.length = 0, segment.colour = .ink$muted, segment.size = 0.25,
      # A white halo keeps a label legible where it crosses a line or a point.
      bg.colour = "white", bg.r = 0.12, max.overlaps = Inf, seed = 42
    )
    # Labels stay above the lowest line any of them belongs to: suggestive
    # labels may sit between the two lines, beside their own peaks, rather
    # than being packed in above the significance line among the towers.
    above <- (if (any(lab$logp < sig)) sugg else sig) + 0.03 * ymax
    p <- p + if (label_angle == 0) {
      do.call(ggrepel::geom_text_repel, c(common, list(
        direction = "both", nudge_y = s * 0.04 * ymax,
        ylim = if (down) c(NA, -above) else c(above, NA),
        box.padding = 0.4, point.padding = 0.15, force = 3)))
    } else {
      do.call(ggrepel::geom_text_repel, c(common, list(
        angle = label_angle, hjust = if (down) 1 else 0, direction = "x",
        nudge_y = s * 0.03 * ymax, box.padding = 0.1)))
    }
  }

  # Threshold P-values are named on a right-hand axis, level with their lines,
  # where no data point can cover them.
  lines_at <- c(sig, sugg)
  lines_lab <- parse(text = vapply(10^-lines_at, .p_plotmath, ""))
  p +
    .x_scale(layout, labels = x_axis, top_title = top_title) +
    ggplot2::scale_y_continuous(
      breaks = s * yb$at, labels = yb$labels,
      limits = if (down) c(-ymax, 0) else c(0, ymax),
      expand = ggplot2::expansion(mult = if (down) c(0, floor_gap) else c(floor_gap, 0)),
      sec.axis = ggplot2::dup_axis(name = NULL, breaks = s * lines_at, labels = lines_lab)
    ) +
    ggplot2::labs(x = if (x_axis) "Chromosome" else NULL,
                  y = expression(-log[10](italic(P)))) +
    theme_maffect(base_size = base_size) + .chr_axis_text(base_size) +
    ggplot2::theme(axis.text.y.right = ggplot2::element_text(colour = .ink$muted,
                                                             size = base_size * .text_size$note,
                                                             hjust = 0)) +
    if (!x_axis) ggplot2::theme(axis.text.x = ggplot2::element_blank())
}

# The effect track. "magnitude" draws |beta| +/- 95% CI, i.e. the effect per
# trait-increasing allele (for odds ratios, the OR of the risk allele): the
# sign of a single variant's effect depends only on which allele the GWAS
# reported, so at a glance the size is what matters. "signed" keeps the sign
# (oriented to the minor allele where an allele frequency is available).
#
# "mirror" (miami_plot()) draws two GWAS in one panel: |beta| of the rows with
# `side` = 1 rises from the zero line, |beta| of those with `side` = -1 hangs
# below it, on one shared axis.
.effect_panel <- function(eff, cloud, layout, cs, mode, scale, log_or, allele,
                          down, base_size, point_size) {
  mirror <- mode == "mirror"
  if (mode == "magnitude" || mirror) {
    side <- if (mirror) eff$side else 1
    m <- abs(eff$beta)
    eff$y <- side * m
    eff$lo <- side * pmax(m - 1.96 * eff$se, 0)
    eff$hi <- side * (m + 1.96 * eff$se)
  } else {
    eff$y <- eff$beta
    eff$lo <- eff$beta - 1.96 * eff$se
    eff$hi <- eff$beta + 1.96 * eff$se
  }
  if (is.null(cloud)) {
    cloud <- data.frame(x = numeric(0), beta = numeric(0), col = eff$col[0],
                        stack = numeric(0), side = numeric(0),
                        tier = factor(character(0), levels = c("base", "suggestive", "gws")))
  }
  cloud$y <- if (mirror) {
    cloud$side * abs(cloud$beta)
  } else if (mode == "magnitude") abs(cloud$beta) else cloud$beta
  # Room for every estimate, but don't let one very wide interval squash the
  # rest; intervals beyond the limits are clipped at the edge.
  top <- max(min(max(abs(c(eff$lo, eff$hi))), 1.6 * max(abs(eff$y))),
             abs(cloud$y))
  lim <- if (mode == "magnitude") c(0, top) else c(-top, top)
  ax <- .effect_axis(abs(eff$y[eff$lead & eff$strong]), top, scale,
                     max_estimate = max(abs(c(eff$y, cloud$y))))
  breaks <- .track_breaks(ax, top, signed = mode != "magnitude", log_or = log_or)
  breaks <- breaks[breaks >= lim[1]]
  # Hanging from the axis, the track's zero is the Manhattan's zero: label once.
  if (down) breaks <- breaks[breaks != 0]
  pts <- point_size * base_size / 10
  eff <- eff[order(eff$stack), ]
  strong <- eff[eff$strong, ]
  weak <- eff[!eff$strong, ]
  flagged <- strong[!is.null(strong$flag) & strong$flag, ]

  ylab <- if (mode != "signed") {
    if (log_or) "OR (risk allele)" else quote(group("|", beta, "|"))
  } else {
    allele_txt <- if (identical(allele, "minor")) "(minor allele)" else "(effect allele)"
    if (log_or) paste("OR", allele_txt) else bquote(beta ~ .(allele_txt))
  }

  cloud <- cloud[order(cloud$tier, cloud$stack), ]
  pts_small <- 0.95 * pts

  p <- ggplot2::ggplot(eff) + .band_layer(layout)
  # The same broken-axis mark as on the P-value axis.
  p <- p + .effect_break(ax, top, signed = mode != "magnitude", axis = "y")
  p +
    ggplot2::geom_hline(yintercept = 0, colour = .ink$axis, linewidth = 0.4) +
    ggplot2::geom_point(data = cloud[cloud$tier == "base", ],
                        ggplot2::aes(x = .data$x, y = .data$y),
                        colour = .ink$noise, size = 0.6 * pts, stroke = 0) +
    ggplot2::geom_point(data = cloud[cloud$tier == "suggestive", ],
                        ggplot2::aes(x = .data$x, y = .data$y, colour = .data$col),
                        shape = 1, stroke = 0.35 * pts, size = pts_small) +
    ggplot2::geom_point(data = cloud[cloud$tier == "gws", ],
                        ggplot2::aes(x = .data$x, y = .data$y, colour = .data$col),
                        shape = 16, size = pts_small, alpha = 0.8) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$x, xend = .data$x, y = 0, yend = .data$y),
                          colour = "#cfcec7", linewidth = 0.35) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$x, xend = .data$x, y = .data$lo,
                                       yend = .data$hi, colour = .data$col),
                          linewidth = 0.8, alpha = 0.55, lineend = "round") +
    ggplot2::geom_point(data = weak, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$col),
                        shape = 21, fill = "white", stroke = 0.55 * pts, size = 1.8 * pts) +
    ggplot2::geom_point(data = strong, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$col),
                        shape = 21, colour = "white", stroke = 0.3 * pts, size = 2.1 * pts) +
    # A lead variant whose effect disagrees with its own case/control allele
    # frequencies: a black ring around its point (see find_artifacts()).
    (if (nrow(flagged)) ggplot2::geom_point(
      data = flagged, ggplot2::aes(x = .data$x, y = .data$y),
      shape = 1, colour = .ink$primary, stroke = 0.5 * pts, size = 3.7 * pts)) +
    (if (nrow(strong)) cs$fill(legend = FALSE)) +   # no fill data, no fill scale
    cs$colour() +
    .x_scale(layout) +
    ggplot2::scale_y_continuous(
      transform = .track_transform(ax$trans, down),
      limits = lim, oob = scales::oob_squish, breaks = breaks,
      # Both halves of a Miami plot are magnitudes: label them as such.
      labels = function(b) .as_plotmath(.effect_labels(if (mirror) abs(b) else b, log_or)),
      expand = ggplot2::expansion(mult = c(0.06, 0.06))
    ) +
    ggplot2::labs(x = "Chromosome", y = ylab) +
    theme_maffect(base_size = base_size) + .chr_axis_text(base_size) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

# The effect axis transform, optionally flipped so the track hangs down from
# its zero line.
.track_transform <- function(trans, down = FALSE) {
  if (!down) return(trans)
  scales::new_transform(
    paste0(trans$name, "-down"),
    transform = function(x) -trans$transform(x),
    inverse = function(x) trans$inverse(-x)
  )
}
