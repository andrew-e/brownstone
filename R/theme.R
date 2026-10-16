# Colours and theme.
#
# Ink and surface colours are one step off white so that the data, not the
# chrome, is the loudest thing on the page.
#
# Minor allele frequency classes default to three distinct hues (blue, aqua,
# orange) rather than one ordered ramp: the job of the encoding is to make rare
# signals jump out, and orange does that where dark blue among light blue
# doesn't. The three colours were checked for colour-vision-deficiency
# separation across all pairs (worst deutan delta-E 9.2, OKLab x100). The
# ordinal ramp (one hue, light = common, dark = rare) is kept as an option and
# is used automatically for four classes, where no three-hue set extends safely.

.ink <- list(
  primary = "#0b0b0b",
  secondary = "#52514e",
  muted = "#898781",
  grid = "#e9e8e3",
  axis = "#c3c2b7",
  band = "#f5f5f2",
  noise = "#c2c1ba",
  noise_alt = "#dad9d3",
  # Significant variants with no frequency: a light, recessive grey (drawn
  # underneath the rest). A mid grey would be confused with aqua under
  # deuteranopia (delta-E 5.1); this one keeps every pair at delta-E >= 9.2.
  unknown = "#bdbcb5"
)

#' Colours for minor allele frequency classes
#'
#' `"contrast"` (the default) gives each class its own hue: blue for common
#' (>= 5%), aqua for low-frequency (1-5%), orange for rare (0.1-1%) and deep
#' berry for ultra-rare (< 0.1%) variants, so rare signals stand out. All
#' pairs, and the light grey used for unknown frequencies, stay distinguishable
#' under protanopia, deuteranopia and tritanopia (worst pair delta-E 9.2, OKLab
#' x100). `"ordinal"` uses one hue, light for common and dark for rare, which
#' reads as an ordered scale but makes rare variants less conspicuous.
#'
#' @param n Number of classes (1 to 4).
#' @param palette `"contrast"` or `"ordinal"`.
#' @return A character vector of hex colours, most common class first.
#' @export
#' @examples
#' maf_colours()
#' maf_colours(palette = "ordinal")
maf_colours <- function(n = 4, palette = c("contrast", "ordinal")) {
  palette <- match.arg(palette)
  if (palette == "ordinal") {
    return(switch(as.character(n),
                  `1` = "#2a78d6",
                  `2` = c("#6da7ec", "#104281"),
                  `3` = c("#6da7ec", "#2a78d6", "#104281"),
                  `4` = c("#86b6ef", "#3987e5", "#1c5cab", "#0d366b"),
                  stop("`maf_colours()` provides 1 to 4 colours.", call. = FALSE)))
  }
  switch(as.character(n),
         `1` = "#2a78d6",
         `2` = c("#2a78d6", "#eb6834"),
         `3` = c("#2a78d6", "#1baf7a", "#eb6834"),
         `4` = c("#2a78d6", "#1baf7a", "#eb6834", "#8f2a5e"),
         stop("`maf_colours()` provides 1 to 4 colours.", call. = FALSE))
}

#' A clean theme for genome-wide plots
#'
#' Minimal chrome: no panel border, no vertical grid, faint horizontal grid,
#' secondary-ink text.
#'
#' @param base_size Base font size, in points.
#' @param base_family Base font family.
#' @return A ggplot2 theme.
#' @export
theme_brownstone <- function(base_size = 10, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = .ink$grid, linewidth = 0.3),
      axis.text = ggplot2::element_text(colour = .ink$secondary),
      axis.title = ggplot2::element_text(colour = .ink$secondary),
      axis.ticks = ggplot2::element_blank(),
      legend.position = "top",
      legend.justification = "left",
      legend.title = ggplot2::element_text(colour = .ink$secondary, size = base_size * 0.9),
      legend.text = ggplot2::element_text(colour = .ink$secondary, size = base_size * 0.9),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.spacing = ggplot2::unit(base_size * 0.4, "pt"),
      plot.title = ggplot2::element_text(colour = .ink$primary, size = base_size * 1.3,
                                         margin = ggplot2::margin(b = base_size * 0.4)),
      plot.subtitle = ggplot2::element_text(colour = .ink$secondary, size = base_size * 1.0),
      plot.caption = ggplot2::element_text(colour = .ink$muted, size = base_size * 0.8),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(6, 10, 4, 6)
    )
}

#' Colour ramps for a continuous minor allele frequency scale
#'
#' Used by `brownstone(maf_scale = "continuous")`. Every ramp changes
#' lightness monotonically, so it stays readable with any form of colour
#' blindness, and puts rare variants at the dark, salient end.
#'
#' @param palette `"viridis"` (the default; purple for rare, green for common),
#'   `"mako"` (deep blue to teal) or `"blue"` (one hue, light to dark).
#' @return A character vector of hex colours, most common first.
#' @export
#' @examples
#' maf_ramp()
maf_ramp <- function(palette = c("viridis", "mako", "blue")) {
  palette <- match.arg(palette)
  switch(palette,
         viridis = c("#7AD151", "#35B779", "#1F988B", "#2A788E", "#39568C", "#472F7D", "#440154"),
         mako = c("#42B9AD", "#359CAA", "#347FA4", "#38639D", "#404589", "#3B2E5D", "#2A1B32"),
         blue = c("#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"))
}
