#' Regional Brownstone plot, with linkage disequilibrium
#'
#' @description
#' Zooms in on one region of the genome, in the spirit of LocusZoom: every
#' variant's -log10(P) by position, coloured by its linkage disequilibrium
#' (LD) with the lead variant, above an effect track and a gene track. Unlike a
#' classic regional plot it keeps what a [brownstone()] plot adds:
#'
#' * **LD with the lead is the colour**: r² in five bins (`ld_breaks`), from
#'   light (no LD) to dark (r² >= 0.8), fetched from LDlink (see [ld_proxy()])
#'   or taken from a table of your own. Variants with no LD information are
#'   drawn hollow. A variant missing from LDlink is assigned to the lowest bin
#'   only when its allele frequency alone caps its r² with the lead below the
#'   first break.
#' * **Minor allele frequency is the shape** of every variant above the
#'   suggestive line: circle (>= 5%), square (1-5%), triangle (0.1-1%) or
#'   star (< 0.1%), drawn with 'ggstar'; a plus marks an unknown frequency.
#'   Each departs further from the plain circle, so the rarer the variant the
#'   more it stands out. A rare variant always has a low
#'   r² with a common lead, so the shape is what shows that a light point is a
#'   rare signal.
#' * **The effect track shows |beta| of every variant above the suggestive
#'   line**, coloured and shaped the same way. The lead, and the strongest
#'   variant of each other frequency class when its signal is too strong to be
#'   an LD shadow of the lead, are drawn with their 95% confidence interval. A
#'   light point well above the lead's |beta| (dotted line) is a second signal.
#' * **The gene track** shows protein-coding genes (Ensembl), the gene at the
#'   lead variant in dark ink.
#'
#' LD is fetched from LDlink by default, which needs a token in `LDLINK_TOKEN`
#' (see [ld_proxy()]); without one the call stops. Without LD
#' (`use_ld_info = FALSE`, or LDlink unreachable) the plot is a zoomed-in
#' Brownstone plot: hits coloured by frequency class.
#'
#' The region is given by one of `snp`, `gene`, or `chr`, `start` and `end`;
#' with none of them, it is centred on the most significant variant. Every
#' call standardises `data`: for several plots of a large GWAS, standardise it
#' once with [standardise_sumstats()] and pass the result.
#'
#' @inheritParams brownstone
#' @param snp The lead variant, by ID (`"rs7412"`, or the data's own IDs) or
#'   position (`"19:44908822"`). The region is `flank` base pairs either side.
#' @param gene A gene symbol: the region is the gene, plus `flank` either side.
#' @param chr,start,end The region by position.
#' @param flank Base pairs either side of `snp` or `gene`, or of the most
#'   significant variant.
#' @param use_ld_info Colour by LD with the lead? `TRUE` (the default) needs LD
#'   from `ld`; `FALSE` draws without LD, ignoring `ld`.
#' @param ld Where r² with the lead comes from: `"ldlink"` (the default; needs a
#'   token, see [ld_proxy()]), a data frame (`snp` or `chr` + `pos`, and `r2`;
#'   or PLINK `--r2` output with `SNP_A`, `CHR_B`, `BP_B`, `SNP_B`, `R2`, of
#'   which the rows for the lead are used), or the name of a column of `data`.
#' @param pop,token The LDlink reference population and access token.
#' @param ld_breaks Boundaries of the r² bins.
#' @param effect_scale How the effect axis reacts to the data, as in
#'   [brownstone()].
#' @param maf_breaks Boundaries of the minor allele frequency classes (at most
#'   four classes).
#' @param genes `TRUE` for the bundled gene table of `build`, a data frame of
#'   your own (`chr`, `start`, `end`, `gene`, optionally `strand`), or `FALSE`
#'   for no gene track.
#' @param max_gene_rows At most this many rows of genes; genes that don't fit
#'   are counted at the right.
#'
#' @return A patchwork object.
#' @export
#' @examples
#' gwas <- simulate_gwas(n_variants = 2e5)
#' # LD as simulated, for the locus at HMGCR:
#' causal <- attr(gwas, "causal")
#' hmgcr <- causal[causal$gene == "HMGCR", ]
#' region_plot(gwas, snp = paste0(hmgcr$chr, ":", hmgcr$pos), ld = attr(gwas, "ld"))
#' \donttest{
#' # LD from LDlink (needs LDLINK_TOKEN):
#' if (nzchar(Sys.getenv("LDLINK_TOKEN"))) region_plot(gwas, gene = "HMGCR")
#' }
#' # No LD:
#' region_plot(gwas, gene = "HMGCR", use_ld_info = FALSE)
region_plot <- function(data, ..., snp = NULL, gene = NULL, chr = NULL, start = NULL,
                        end = NULL, flank = 2.5e5, build = c("GRCh38", "GRCh37"),
                        use_ld_info = TRUE, ld = "ldlink", pop = "EUR",
                        token = Sys.getenv("LDLINK_TOKEN"),
                        threshold = 5e-8, suggestive = 1e-5,
                        ld_breaks = c(0.2, 0.4, 0.6, 0.8),
                        effect_scale = c("auto", "sqrt", "linear"),
                        maf_breaks = c(0.001, 0.01, 0.05), genes = TRUE, max_gene_rows = 6,
                        base_size = 10, point_size = 1) {
  build <- match.arg(build)
  effect_scale <- match.arg(effect_scale)
  if ((!is.null(snp)) + (!is.null(gene)) + (!is.null(chr) || !is.null(start) || !is.null(end)) > 1) {
    stop("Give the region as one of `snp`, `gene`, or `chr`, `start` and `end`.", call. = FALSE)
  }
  if (use_ld_info && identical(ld, "ldlink") && !nzchar(token)) {
    stop("LDLINK_TOKEN is not set. If you want LD information, get a token here: ",
         "https://ldlink.nih.gov/apiaccess and add LDLINK_TOKEN=<token> to your .Renviron; ",
         "else set `use_ld_info = FALSE` to continue without.", call. = FALSE)
  }
  ld_col <- if (use_ld_info && is.character(ld) && length(ld) == 1L && ld != "ldlink") ld
  d <- .standardise(data, ..., keep = ld_col, ids = FALSE)
  sig <- -log10(threshold)
  if (!is.null(suggestive) && suggestive <= threshold) {
    stop("`suggestive` must be a larger P-value than `threshold`.", call. = FALSE)
  }
  sugg <- if (!is.null(suggestive)) -log10(suggestive)
  floor_logp <- sugg %||% sig
  gene_table <- if (is.data.frame(genes)) genes

  reg <- .region_window(d, snp, gene, chr, start, end, flank, build, gene_table)
  rows <- which(d$chr == reg$chr & d$pos >= reg$start & d$pos <= reg$end)
  if (!length(rows)) {
    stop("No variants in ", reg$chr, ":", reg$start, "-", reg$end, ".", call. = FALSE)
  }
  r <- d[rows, ]
  li <- if (!is.null(reg$lead)) match(reg$lead, rows) else which.max(r$logp)
  lead_id <- .variant_id(r, li)
  has_maf <- isTRUE(attr(d, "has_maf"))

  # LD with the lead, binned; NULL without LD.
  r2 <- if (use_ld_info) .region_ld(ld, ld_col, r, li, lead_id, reg, build, pop, token)
  ldc <- if (!is.null(r2)) .ld_classes(r2, r$maf, r$maf[li], ld_breaks, has_maf)
  if (!is.null(ldc) && has_maf && isTRUE(r$maf[li] < 0.01)) {
    message("The lead variant is rare (MAF ", signif(100 * r$maf[li], 2), "%): r\u00b2 from ",
            "1000 Genomes is unreliable at that frequency.")
  }

  r$class <- if (has_maf) .maf_classes(r$maf, maf_breaks) else factor(rep("all", nrow(r)))
  if (anyNA(r$class)) {
    r$class <- factor(r$class, levels = c(levels(r$class), "Unknown"))
    r$class[is.na(r$class)] <- "Unknown"
  }
  r$y_logp <- r$logp
  r$hit <- r$logp >= floor_logp
  r$lead <- seq_len(nrow(r)) == li
  pts <- point_size * base_size / 10

  has_effect <- isTRUE(attr(d, "has_effect")) && is.finite(r$beta[li]) && is.finite(r$se[li])
  main <- .region_main(r, ldc, reg, sig, sugg, floor_logp, lead_id, has_maf, maf_breaks,
                       base_size, pts, floor_gap = if (has_effect) 0.03 else 0)
  if (!has_effect) message("No effect size and standard error for the lead variant; drawing without the effect track.")
  track <- if (has_effect) {
    .region_track(r, ldc, reg, r2, effect_scale, log_or = identical(attr(d, "effect_scale"), "log_or"),
                  allele = attr(d, "effect_allele"), has_maf, maf_breaks, base_size, pts)
  }
  gene_panel <- if (!isFALSE(genes)) {
    .region_genes(reg, r$pos[li], build, gene_table, max_gene_rows, base_size)
  }
  panels <- Filter(Negate(is.null), list(main, track, gene_panel$plot))
  heights <- c(1, if (!is.null(track)) 0.35, gene_panel$height)
  # Only the bottom panel carries the position axis; a bold ground line and a
  # little air separate the -log10(P) panel from the effect track.
  for (k in seq_along(panels)) {
    last <- k == length(panels)
    top <- if (k == 1) 6 else if (k == 2 && !is.null(track)) 8 else 2
    panels[[k]] <- panels[[k]] + .region_x(reg, axis = last) +
      (if (k == 1 && !is.null(track)) .ground_rule()) +
      ggplot2::theme(plot.margin = ggplot2::margin(top, 10, if (last) 4 else 0, 6))
  }
  # Titles added with plot_annotation() take the theme's (large, centred) style.
  patchwork::wrap_plots(panels, ncol = 1, heights = heights) +
    patchwork::plot_annotation(theme = theme_brownstone(base_size = base_size))
}

# The region and, if given by `snp`, its lead variant (a row of `d`).
.region_window <- function(d, snp, gene, chr, start, end, flank, build, gene_table) {
  if (!is.null(snp)) {
    snp <- as.character(snp)
    m <- regmatches(snp, regexec("^(chr)?([0-9]{1,2}|X|Y|XY|MT)[:_]([0-9]+)$", snp,
                                 ignore.case = TRUE))[[1]]
    i <- if (length(m)) {
      at <- which(d$chr == .norm_chr(m[3]) & d$pos == as.numeric(m[4]))
      at[which.max(d$logp[at])]
    } else if ("snp" %in% names(attr(d, "columns"))) {
      match(snp, d$snp)
    }
    if (!length(i) || is.na(i)) stop("Variant ", snp, " is not in the data.", call. = FALSE)
    return(list(chr = d$chr[i], start = d$pos[i] - flank, end = d$pos[i] + flank, lead = i))
  }
  if (!is.null(gene)) {
    g <- gene_table %||% .gene_table(build)
    g <- g[toupper(g$gene) == toupper(gene), , drop = FALSE]
    if (!nrow(g)) stop("Gene ", gene, " is not in the ", build, " gene table.", call. = FALSE)
    return(list(chr = .norm_chr(g$chr[1]), start = g$start[1] - flank, end = g$end[1] + flank))
  }
  if (!is.null(chr) || !is.null(start) || !is.null(end)) {
    if (is.null(chr) || is.null(start) || is.null(end) || start >= end) {
      stop("Give `chr`, `start` and `end` together, with `start` < `end`.", call. = FALSE)
    }
    return(list(chr = .norm_chr(chr), start = start, end = end))
  }
  i <- which.max(d$logp)
  list(chr = d$chr[i], start = d$pos[i] - flank, end = d$pos[i] + flank, lead = i)
}

# r² of every variant of the region `r` with its lead (row `li`), or NULL.
.region_ld <- function(ld, ld_col, r, li, lead_id, reg, build, pop, token) {
  if (!is.null(ld_col)) {
    r2 <- as.numeric(r[[ld_col]])
    r2[li] <- 1
    return(r2)
  }
  source <- "your table"
  by_id <- FALSE
  if (is.data.frame(ld)) {
    tab <- .ld_table(ld, lead_id, r$chr[li], r$pos[li])
    if (is.null(tab)) return(NULL)
  } else if (identical(ld, "ldlink")) {
    by_id <- grepl("^rs[0-9]+$", lead_id, ignore.case = TRUE)
    var <- if (by_id) lead_id else paste0(r$chr[li], ":", r$pos[li])
    window <- min(1e6, max(r$pos[li] - reg$start, reg$end - r$pos[li]))
    tab <- tryCatch(ld_proxy(var, pop = pop, build = build, window = window, token = token),
                    error = function(e) {
                      message(conditionMessage(e), " Drawing without LD.")
                      NULL
                    })
    if (is.null(tab)) return(NULL)
    source <- paste0("LDlink, ", paste(pop, collapse = "+"))
    # The query's own row says where LDlink thinks the lead is: a different
    # position means the builds differ, and positions can't be trusted.
    q <- tab[which(tab$distance == 0 & tolower(tab$snp) == tolower(lead_id))[1], ]
    if (by_id && nrow(q) && !is.na(q$pos) && q$pos != r$pos[li]) {
      warning("LDlink places ", lead_id, " at ", q$chr, ":", q$pos, ", not ", r$chr[li], ":",
              r$pos[li], " as in the data: is `build` right? Matching LD by rsID only.",
              call. = FALSE)
      tab$chr <- NULL
      tab$pos <- NULL
    }
  } else {
    stop("`ld` must be \"ldlink\", a data frame or a column name; ",
         "for no LD, set `use_ld_info = FALSE`.", call. = FALSE)
  }
  i <- .match_variants(tab, r, unique_pos = TRUE)
  ok <- !is.na(i)
  if (!any(ok)) {
    message("None of the LD (", source, ") matches a variant in the region; drawing without LD.")
    return(NULL)
  }
  if (identical(ld, "ldlink") && !by_id && sum(ok) < 0.1 * nrow(tab)) {
    warning("Only ", sum(ok), " of ", nrow(tab), " LDlink variants are at positions in the data: ",
            "is `build` right?", call. = FALSE)
  }
  r2 <- rep(NA_real_, nrow(r))
  r2[i[ok]] <- tab$r2[ok]
  r2[li] <- 1
  message("LD with ", lead_id, " (", source, "): r\u00b2 for ", sum(!is.na(r2)), " of ",
          nrow(r), " variants in the region.")
  r2
}

# A table of LD in standard columns (snp, chr, pos, r2), from one with columns
# named loosely or from PLINK --r2 output (rows for the lead only).
.ld_table <- function(ld, lead_id, lead_chr, lead_pos) {
  nm <- tolower(names(ld))
  col <- function(...) {
    k <- match(c(...), nm)
    k <- k[!is.na(k)]
    if (length(k)) ld[[k[1]]]
  }
  if (!is.null(col("snp_a"))) {
    keep <- as.character(col("snp_a")) %in% lead_id
    if (!is.null(col("bp_a")) && !is.null(col("chr_a"))) {
      keep <- keep | (.norm_chr(col("chr_a")) == lead_chr & as.numeric(col("bp_a")) == lead_pos)
    }
    keep[is.na(keep)] <- FALSE
    if (!any(keep)) {
      message("The LD table has no rows for ", lead_id, ", so no LD.")
      return(NULL)
    }
    out <- list(snp = col("snp_b")[keep], chr = col("chr_b")[keep], pos = col("bp_b")[keep],
                r2 = col("r2")[keep])
  } else {
    out <- list(snp = col("snp", "rsid", "id", "variant"), chr = col("chr", "chrom", "chromosome"),
                pos = col("pos", "bp", "position"), r2 = col("r2", "rsq", "r_squared", "r.squared"))
  }
  if (is.null(out$r2)) stop("An LD table needs an `r2` column.", call. = FALSE)
  if (is.null(out$snp) && (is.null(out$chr) || is.null(out$pos))) {
    stop("An LD table needs a `snp` column, or `chr` and `pos` columns.", call. = FALSE)
  }
  out <- as.data.frame(Filter(Negate(is.null), out), stringsAsFactors = FALSE)
  out$r2 <- as.numeric(out$r2)
  out
}

# r² bins: the lowest also takes variants whose allele frequencies alone cap
# their r² with the lead below the first break; other variants with no r² are
# "No LD data" (drawn hollow).
.ld_classes <- function(r2, maf, lead_maf, breaks, has_maf) {
  n <- length(breaks)
  num <- function(b) format(b, drop0trailing = TRUE, trim = TRUE)
  labs <- c(paste0("< ", num(breaks[1])), paste0(num(breaks[-n]), "-", num(breaks[-1])),
            paste0(">= ", num(breaks[n])))
  cls <- cut(r2, c(-Inf, breaks, Inf), labels = labs, right = FALSE)
  if (has_maf && is.finite(lead_maf)) {
    capped <- is.na(r2) & is.finite(maf) & .r_max(maf, lead_maf)^2 < breaks[1]
    cls[capped] <- labs[1]
  }
  cls <- factor(cls, levels = c(labs, "No LD data"))
  cls[is.na(cls)] <- "No LD data"
  cls
}

# Colours for the r² bins, light (no LD) to dark: one hue (violet) with
# lightness falling monotonically, so the order reads under any colour vision.
# Violet is the hue the frequency-class colours don't use: with normal colour
# vision every bin is at least delta-E 13 from all of them, so a colour in a
# regional plot is never mistaken for a frequency class in another plot. (No
# five-step ramp can avoid all four under every colour-vision deficiency.)
.ld_colours <- function(n = 5) {
  ramp <- c("#c1a8f0", "#9974d4", "#7044ac", "#471d77", "#1e043a")
  if (n == length(ramp)) return(ramp)
  scales::pal_gradient_n(ramp)(seq(0, 1, length.out = n))
}

# Shapes (ggstar) for the minor allele frequency classes, most common first:
# circle, square, triangle, star. Each departs further from the plain circle
# of a common variant, so neighbouring classes look alike and the rarest
# stand out most; unknown frequency is a plus, unlike any of them.
.maf_shapes <- function(levels) {
  shapes <- c(15, 13, 11, 1)   # ggstar: circle, square, regular triangle, pentagram
  out <- stats::setNames(rep(15, length(levels)), levels)
  known <- setdiff(levels, "Unknown")
  if (length(known) > length(shapes)) {
    stop("A regional plot shows at most four frequency classes (`maf_breaks` of length 3).",
         call. = FALSE)
  }
  out[known] <- shapes[seq_along(known)]
  if ("Unknown" %in% levels) out["Unknown"] <- 29   # ggstar: plus filled
  out
}

# Size multipliers that make the shapes look equally large: at the same size
# a triangle or a star reads smaller than a circle, a square larger.
.shape_size <- function(shape) {
  c(`15` = 1, `13` = 0.92, `11` = 1.25, `1` = 1.15, `29` = 1.05)[as.character(shape)]
}

# Fill (r² bin) and shape (frequency class) scales shared by both panels.
.region_scales <- function(ldc, r, legend, pts = 1) {
  lev <- levels(ldc)
  fills <- stats::setNames(c(.ld_colours(length(lev) - 1L), "white"), lev)
  shapes <- .maf_shapes(levels(r$class))
  present <- levels(r$class)[levels(r$class) %in% r$class[r$hit | r$lead]]
  list(
    ggplot2::scale_fill_manual(
      values = fills, limits = lev, drop = FALSE, labels = .maf_legend_labels(lev),
      # Just r²: the lead is labelled, and a one-row legend then fits a
      # single-column figure.
      name = quote(italic(r)^2),
      guide = if (legend) ggplot2::guide_legend(order = 1, nrow = 1, override.aes = list(
        shape = 21, starshape = 15, size = 3 * pts, starstroke = 0.3,
        colour = c(rep(.ink$secondary, length(lev) - 1L), .ink$muted)))
      else "none"),
    ggstar::scale_starshape_manual(
      values = shapes, limits = present,
      labels = .maf_legend_labels(present), name = "Minor allele frequency",
      guide = if (legend && length(present)) {
        ggplot2::guide_legend(order = 2, override.aes = list(
          fill = "white", colour = .ink$secondary, starstroke = 0.5,
          size = 2.8 * pts * .shape_size(shapes[present])))
      } else "none")
  )
}

.region_x <- function(reg, axis) {
  list(
    ggplot2::scale_x_continuous(limits = c(reg$start, reg$end), expand = ggplot2::expansion(0),
                                labels = function(b) format(b / 1e6, trim = TRUE)),
    ggplot2::labs(x = if (axis) paste0("Position on chromosome ", reg$chr, " (Mb)")),
    if (!axis) ggplot2::theme(axis.text.x = ggplot2::element_blank())
  )
}

# Variants drawn in LD colours: the background as small circles, hits in their
# frequency-class shape with a dark edge (hollow where there is no LD data),
# rarest last, and the lead large with a black edge.
.ld_points <- function(r, ldc, y, pts) {
  r$y <- y
  r$ld <- ldc
  r$edge <- ifelse(r$ld == "No LD data", .ink$muted, .ink$secondary)
  bg <- r[!r$hit & !r$lead, ]
  hits <- r[r$hit & !r$lead, ]
  hits <- hits[order(ifelse(hits$class %in% "Unknown", 0L, as.integer(hits$class))), ]
  lead <- r[r$lead, ]
  shapes <- .maf_shapes(levels(r$class))
  list(
    ggplot2::geom_point(data = bg, ggplot2::aes(x = .data$pos, y = .data$y, fill = .data$ld),
                        shape = 21,
                        colour = ifelse(bg$ld == "No LD data", .ink$muted, "transparent"),
                        stroke = 0.25 * pts, size = 0.9 * pts),
    ggstar::geom_star(data = hits, ggplot2::aes(x = .data$pos, y = .data$y, fill = .data$ld,
                                                starshape = .data$class),
                      colour = hits$edge, starstroke = 0.3 * pts,
                      size = 2 * pts * .shape_size(shapes[as.character(hits$class)])),
    ggstar::geom_star(data = lead, ggplot2::aes(x = .data$pos, y = .data$y, fill = .data$ld,
                                                starshape = .data$class),
                      colour = .ink$primary, starstroke = 0.6 * pts,
                      size = 3.6 * pts * .shape_size(shapes[as.character(lead$class)]))
  )
}

# Variants drawn without LD, as in brownstone(): a grey floor and hits in
# their frequency-class colour (solid, or rings for suggestive peaks).
.maf_points <- function(r, y, sig, sugg, floor_logp, has_maf, maf_breaks, pts, legend) {
  r$y <- y
  r$x <- r$pos
  cp <- .class_palette(r, has_maf, maf_breaks, NULL, floor_logp)
  r$class <- cp$class
  cs <- .colour_spec(r, has_maf, "binned", cp$pal, NULL, floor_logp)
  r$col <- cs$values
  r$stack <- cs$stack
  hits <- r[r$hit, ]
  hits$tier <- factor(ifelse(hits$logp >= sig, "gws", "suggestive"), levels = c("gws", "suggestive"))
  list(
    ggplot2::geom_point(data = r[!r$hit, ], ggplot2::aes(x = .data$x, y = .data$y),
                        colour = .ink$noise, size = 0.6 * pts, stroke = 0),
    .hit_layers(hits, cs, legend = legend && has_maf, sig, sugg, pts)
  )
}

.region_main <- function(r, ldc, reg, sig, sugg, floor_logp, lead_id, has_maf, maf_breaks,
                         base_size, pts, floor_gap = 0) {
  max_logp <- max(r$logp)
  comp <- .y_compressor(max_logp, sig)
  y <- comp$f(r$logp)
  ymax <- max(comp$f(max_logp) * 1.12, sig * 1.05)
  p <- ggplot2::ggplot()
  if (!is.null(comp$at)) p <- p + .axis_break(comp$at, ymax)
  p <- p + .p_lines(sig, sugg) +
    if (is.null(ldc)) {
      .maf_points(r, y, sig, sugg, floor_logp, has_maf, maf_breaks, pts, legend = TRUE)
    } else {
      c(.ld_points(r, ldc, y, pts), .region_scales(ldc, r, legend = TRUE, pts = pts))
    }
  li <- which(r$lead)
  lab <- data.frame(x = r$pos[li], y = y[li], label = lead_id, fontface = "bold",
                    colour = .ink$primary)
  near <- data.frame(x = r$pos[r$hit], y = y[r$hit])
  p + .repel_labels(lab, near, base_size, nudge_y = 0.06 * ymax) +
    .p_axis(comp, max_logp, ymax, sig, sugg, below = floor_gap) +
    ggplot2::labs(y = expression(-log[10](italic(P)))) +
    theme_brownstone(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(colour = .ink$grid, linewidth = 0.3),
      axis.text.y.right = ggplot2::element_text(colour = .ink$muted, size = base_size * .text_size$note,
                                                hjust = 0),
      legend.box = "vertical", legend.box.just = "left"
    )
}

# The effect track: |beta| of every hit. The lead, and the strongest variant
# of each other frequency class when it is too significant to be an LD shadow
# of the lead (its |z| beyond r * |z_lead| + 2, with r measured where known
# and otherwise the most that allele frequencies allow), get their 95% CI.
.region_track <- function(r, ldc, reg, r2, scale, log_or, allele, has_maf, maf_breaks,
                          base_size, pts) {
  r$y <- abs(r$beta)
  r$lo <- pmax(r$y - 1.96 * r$se, 0)
  r$hi <- r$y + 1.96 * r$se
  on <- (r$hit | r$lead) & is.finite(r$y) & is.finite(r$se)
  li <- which(r$lead)
  ci <- li
  z <- .z_from_logp(r$logp)
  for (cl in setdiff(levels(r$class), c(as.character(r$class[li]), "Unknown"))) {
    k <- which(on & r$class == cl)
    if (!length(k)) next
    k <- k[which.max(r$logp[k])]
    bound <- if (!is.null(r2) && is.finite(r2[k])) sqrt(r2[k]) else .r_max(r$maf[k], r$maf[li])
    if (isTRUE(z[k] > bound * z[li] + 2)) ci <- c(ci, k)
  }
  shown <- r[on, ]
  cis <- r[ci, ]
  top <- max(min(max(cis$hi), 1.6 * max(shown$y)), shown$y)
  ax <- .effect_axis(r$y[li], top, scale, max_estimate = max(shown$y))
  breaks <- .track_breaks(ax, top, signed = FALSE, log_or = log_or)
  p <- ggplot2::ggplot() +
    .effect_break(ax, top, signed = FALSE, axis = "y") +
    ggplot2::geom_hline(yintercept = 0, colour = .ink$axis, linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = r$y[li], colour = .ink$muted, linewidth = 0.3, linetype = "13") +
    ggplot2::geom_segment(data = cis, ggplot2::aes(x = .data$pos, xend = .data$pos, y = 0,
                                                   yend = .data$y),
                          colour = "#cfcec7", linewidth = 0.35) +
    ggplot2::geom_segment(data = cis, ggplot2::aes(x = .data$pos, xend = .data$pos, y = .data$lo,
                                                   yend = .data$hi),
                          colour = .ink$secondary, linewidth = 0.8, alpha = 0.45, lineend = "round")
  p <- p + if (is.null(ldc)) {
    shown$x <- shown$pos
    cp <- .class_palette(shown, has_maf, maf_breaks, NULL, -Inf)
    shown$class <- cp$class
    cs <- .colour_spec(shown, has_maf, "binned", cp$pal, NULL, -Inf)
    list(ggplot2::geom_point(data = shown[order(cs$stack), ],
                             ggplot2::aes(x = .data$pos, y = .data$y, fill = .data$class),
                             shape = 21, colour = "white", stroke = 0.3 * pts,
                             size = ifelse(shown$lead[order(cs$stack)], 3.2, 2) * pts),
         cs$fill(legend = FALSE))
  } else {
    c(.ld_points(r[on, ], ldc[on], r$y[on], pts),
      .region_scales(ldc, r, legend = FALSE))
  }
  p +
    ggplot2::scale_y_continuous(
      transform = ax$trans, limits = c(0, top), oob = scales::oob_squish, breaks = breaks,
      labels = function(b) .as_plotmath(.effect_labels(b, log_or)),
      expand = ggplot2::expansion(mult = c(0.06, 0.08))
    ) +
    ggplot2::labs(y = .effect_title(signed = FALSE, log_or, allele)) +
    theme_brownstone(base_size = base_size) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

# The gene track, and its height relative to the association panel.
.region_genes <- function(reg, lead_pos, build, gene_table, max_rows, base_size) {
  g <- .genes_in_window(reg$chr, reg$start, reg$end, build, gene_table)
  g <- .pack_genes(g, reg$start, reg$end, lead_pos, base_size, max_rows)
  shown <- g[!is.na(g$row), ]
  n_rows <- max(c(1, shown$row))
  shown$x0 <- pmax(ifelse(shown$strand == "-", shown$end, shown$start), reg$start)
  shown$x1 <- pmin(ifelse(shown$strand == "-", shown$start, shown$end), reg$end)
  shown$x0 <- pmin(shown$x0, reg$end)
  shown$x1 <- pmax(shown$x1, reg$start)
  at_lead <- shown$start <= lead_pos & shown$end >= lead_pos
  ink <- ifelse(at_lead, .ink$primary, .ink$secondary)
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = shown, ggplot2::aes(x = .data$x0, xend = .data$x1, y = -.data$row,
                                                     yend = -.data$row),
                          colour = ink, linewidth = 0.6,
                          arrow = ggplot2::arrow(length = ggplot2::unit(0.05, "inches"), type = "closed")) +
    # In a figure narrower than the 7 in the spacing was estimated for, a name
    # that would overlap one drawn before it is left out (names are drawn
    # nearest the lead first, so the lead's gene always keeps its name).
    ggplot2::geom_text(data = shown, ggplot2::aes(x = .data$label_x, y = -.data$row - 0.42,
                                                  label = .data$gene, hjust = .data$hjust),
                       colour = ink, fontface = "italic", check_overlap = TRUE,
                       size = base_size * .text_size$note / ggplot2::.pt)
  if (!nrow(g)) {
    p <- p + ggplot2::annotate("text", x = (reg$start + reg$end) / 2, y = -1, colour = .ink$muted,
                               label = "No protein-coding genes",
                               size = base_size * .text_size$note / ggplot2::.pt)
  }
  hidden <- sum(is.na(g$row))
  if (hidden) {
    p <- p + ggplot2::annotate("text", x = reg$end, y = -n_rows - 0.42, hjust = 1, vjust = 0.5,
                               colour = .ink$muted, label = paste0("+", hidden, " genes"),
                               size = base_size * .text_size$note / ggplot2::.pt)
  }
  p <- p +
    ggplot2::scale_y_continuous(limits = c(-n_rows - 0.75, -0.6), expand = ggplot2::expansion(0)) +
    # In a narrow figure, a name near an edge may run into the margin.
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(y = NULL) +
    theme_brownstone(base_size = base_size) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank())
  list(plot = p, height = 0.05 + 0.07 * n_rows)
}
