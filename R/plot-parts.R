# Pieces shared by the scatter-type plots (volcano_plot(), effect_vs_maf()).
# Hits are drawn as in brownstone()'s main panel, so the same variant looks the
# same in every plot of the package.

# Frequency classes and their colours. Significant variants with no frequency
# get their own neutral class rather than silently looking like noise.
.class_palette <- function(d, has_maf, maf_breaks, colours, floor_logp) {
  cls <- if (has_maf) .maf_classes(d$maf, maf_breaks) else factor(rep("all", nrow(d)))
  n <- nlevels(cls)
  pal <- if (has_maf) colours %||% maf_colours(n) else (colours %||% maf_colours(1))[1]
  if (length(pal) < n) stop("`colours` needs ", n, " colours.", call. = FALSE)
  pal <- stats::setNames(pal[seq_len(n)], levels(cls))
  if (has_maf && anyNA(cls[d$logp >= floor_logp])) {
    cls <- factor(cls, levels = c(levels(cls), "Unknown"))
    cls[is.na(cls)] <- "Unknown"
    pal <- c(pal, Unknown = .ink$unknown)
  }
  list(class = cls, pal = pal)
}

# Hits: genome-wide significant ones solid, suggestive peaks as hollow rings,
# both in their colour (`hits$col`, scaled by `cs`). `shape` is mapped to the
# tier only to get a legend entry for it, when there is a suggestive line.
.hit_layers <- function(hits, cs, legend, sig, sugg, pts) {
  if (!nrow(hits)) return(NULL)
  hits <- hits[order(hits$stack), ]   # rarest last, i.e. on top
  strong <- hits[hits$tier == "gws", ]
  weak <- hits[hits$tier == "suggestive", ]
  list(
    ggplot2::geom_point(data = weak, ggplot2::aes(x = .data$x, y = .data$y,
                                                  colour = .data$col, shape = .data$tier),
                        fill = "white", stroke = 0.6 * pts, size = 1.45 * pts),
    ggplot2::geom_point(data = strong, ggplot2::aes(x = .data$x, y = .data$y,
                                                    fill = .data$col, shape = .data$tier),
                        colour = "white", stroke = 0.25 * pts, size = 1.7 * pts),
    # The colour legend comes from the solid points, or from the rings when
    # nothing is genome-wide significant.
    if (nrow(strong)) cs$fill(legend = legend),
    cs$colour(legend = legend && !nrow(strong)),
    ggplot2::scale_shape_manual(
      values = c(gws = 21, suggestive = 21), limits = c("gws", "suggestive"),
      labels = parse(text = paste("italic(P) <",
                                  vapply(10^-c(sig, sugg %||% sig), .p_plotmath, ""))),
      name = NULL,
      guide = if (!is.null(sugg)) {
        ggplot2::guide_legend(order = 2, override.aes = list(
          fill = c(.ink$secondary, "white"), colour = c("white", .ink$secondary),
          stroke = c(0.25, 0.6), size = c(3, 2.6)))
      } else "none"
    )
  )
}

# Significance lines, as in brownstone(): their P-values are named on a
# right-hand axis (see .p_axis()).
.p_lines <- function(sig, sugg) {
  list(
    if (!is.null(sugg)) {
      ggplot2::geom_hline(yintercept = sugg, colour = .ink$muted, linewidth = 0.25,
                          linetype = "13")
    },
    ggplot2::geom_hline(yintercept = sig, colour = .ink$secondary, linewidth = 0.3,
                        linetype = "22")
  )
}

# A -log10(P) axis (see .y_compressor()), with the significance lines' P-values
# on the right; `below` is space under 0, as a fraction of the axis.
.p_axis <- function(comp, max_logp, ymax, sig, sugg, below = 0) {
  yb <- .y_breaks(comp, max_logp)
  lines_at <- c(sig, sugg)
  ggplot2::scale_y_continuous(
    breaks = yb$at, labels = yb$labels, limits = c(0, ymax),
    expand = ggplot2::expansion(mult = c(below, 0)),
    sec.axis = ggplot2::dup_axis(name = NULL, breaks = lines_at,
                                 labels = parse(text = vapply(10^-lines_at, .p_plotmath, "")))
  )
}

# Locus labels (see .locus_labels()) placed with ggrepel, repelling from the
# points in `near` too: a coarse grid of the hits (thinned by the caller, on
# the axes as drawn) is enough for that and keeps ggrepel fast.
.repel_labels <- function(lab, near, base_size, nudge_y = 0, ylim = c(NA, NA)) {
  if (!nrow(lab)) return(NULL)
  rep <- rbind(lab[, c("x", "y", "label", "fontface", "colour")],
               data.frame(x = near$x, y = near$y, label = character(nrow(near)),
                          fontface = rep("plain", nrow(near)),
                          colour = rep(.ink$primary, nrow(near))))
  ggrepel::geom_text_repel(
    data = rep, mapping = ggplot2::aes(x = .data$x, y = .data$y, label = .data$label,
                                       fontface = .data$fontface),
    colour = rep$colour, size = base_size * .text_size$label / ggplot2::.pt,
    min.segment.length = 0, segment.colour = .ink$muted, segment.size = 0.25,
    max.overlaps = Inf, seed = 42, direction = "both", nudge_y = nudge_y,
    ylim = ylim, box.padding = 0.4, point.padding = 0.15, force = 3
  )
}

# A broken-axis mark across an effect axis (see .effect_axis()) at +-`at`, for
# `axis = "x"` or `"y"`: the same white strip edged with hairlines as on the
# -log10(P) axis. Drawn over a dense cloud of points, it shows where the axis
# changes scale.
.effect_break <- function(ax, top, signed, axis = c("y", "x")) {
  axis <- match.arg(axis)
  if (is.null(ax$at)) return(NULL)
  gap <- 0.025 * ax$trans$transform(top)
  edges <- c(ax$at - gap, ax$trans$inverse(ax$at + gap))
  sides <- if (signed) c(1, -1) else 1
  strips <- lapply(sides, function(s) {
    e <- sort(s * edges)
    if (axis == "y") {
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = e[1], ymax = e[2],
                        fill = "white", colour = NA)
    } else {
      ggplot2::annotate("rect", ymin = -Inf, ymax = Inf, xmin = e[1], xmax = e[2],
                        fill = "white", colour = NA)
    }
  })
  at <- as.vector(outer(edges, sides))
  hair <- if (axis == "y") {
    ggplot2::geom_hline(yintercept = at, colour = .ink$axis, linewidth = 0.3)
  } else {
    ggplot2::geom_vline(xintercept = at, colour = .ink$axis, linewidth = 0.3)
  }
  c(strips, list(hair))
}

# A bold rule along the foot of a -log10(P) panel that has an effect track
# beneath it: the ground the towers stand on, so the two panels read apart.
.ground_rule <- function() {
  ggplot2::theme(axis.line.x.bottom = ggplot2::element_line(colour = .ink$secondary,
                                                            linewidth = 1.1))
}

# Text as a plotmath string constant, quotes escaped.
.pm_string <- function(x) sprintf("'%s'", gsub("'", "\\\\'", x))

# Axis title for effects: |beta| or the OR of the risk allele (magnitudes), or
# beta / OR per minor (or effect) allele (signed).
.effect_title <- function(signed, log_or, allele) {
  if (!signed) return(if (log_or) "OR (risk allele)" else quote(group("|", beta, "|")))
  allele_txt <- if (identical(allele, "minor")) "(minor allele)" else "(effect allele)"
  if (log_or) paste("OR", allele_txt) else bquote(beta ~ .(allele_txt))
}
