# Fonts. Every text the package draws names one family, "maffect-sans"
# (or the `maffect.family` option). When the package loads, that name is
# registered with every graphics device that takes font names: devices that
# read installed fonts (ragg, svglite, quartz, windows) get the best modern
# sans already on the system, and pdf() and postscript(), which only know
# their built-in Type 1 fonts, get Helvetica's, so they never fail on the name.
# Cairo devices that don't know the name fall back to their default sans.

.font_alias <- "maffect-sans"

# In order of preference: two widely installed open fonts, then the modern
# sans each operating system ships (macOS, Windows, then Linux).
.font_candidates <- c("Source Sans 3", "Inter", "Avenir Next", "Segoe UI", "Noto Sans",
                      "Open Sans")

# The family every text layer uses.
.family <- function() getOption("maffect.family", .font_alias)

# The first candidate installed on this system, or NULL. systemfonts always
# returns some font, so a match is one whose family is the one asked for.
.find_sans <- function(candidates = .font_candidates) {
  for (f in candidates) {
    info <- tryCatch(systemfonts::font_info(family = f), error = function(e) NULL)
    if (!is.null(info) && nrow(info) && identical(info$family[1], f)) return(f)
  }
  NULL
}

.register_fonts <- function(alias = .font_alias) {
  sans <- .find_sans()
  # ragg, svglite and textshaping find fonts through systemfonts.
  try(systemfonts::register_variant(alias, sans %||% "sans"), silent = TRUE)
  # pdf() and postscript(): Helvetica's metrics under the alias.
  if (is.null(grDevices::pdfFonts(alias)[[1]])) {
    do.call(grDevices::pdfFonts, stats::setNames(grDevices::pdfFonts("Helvetica"), alias))
  }
  if (is.null(grDevices::postscriptFonts(alias)[[1]])) {
    do.call(grDevices::postscriptFonts,
            stats::setNames(grDevices::postscriptFonts("Helvetica"), alias))
  }
  if (is.null(sans)) return(invisible(NULL))
  # The on-screen and bitmap devices of macOS and Windows.
  if (exists("quartzFonts", envir = asNamespace("grDevices"))) {
    styles <- list(c("normal", FALSE), c("bold", FALSE), c("normal", TRUE), c("bold", TRUE))
    ps <- vapply(styles, function(s) {
      systemfonts::font_info(family = sans, weight = s[1], italic = as.logical(s[2]))$name[1]
    }, "")
    try(do.call(grDevices::quartzFonts,
                stats::setNames(list(grDevices::quartzFont(ps)), alias)), silent = TRUE)
  }
  if (.Platform$OS.type == "windows") {
    wf <- get("windowsFonts", envir = asNamespace("grDevices"))
    wfont <- get("windowsFont", envir = asNamespace("grDevices"))
    try(do.call(wf, stats::setNames(list(wfont(sans)), alias)), silent = TRUE)
  }
  invisible(sans)
}

.onLoad <- function(libname, pkgname) {
  try(.register_fonts(), silent = TRUE)
}
