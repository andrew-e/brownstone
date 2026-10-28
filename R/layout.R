# Genome layout, point thinning and axis transforms.

# One row per chromosome present in `d`, in karyotype order, with its offset on
# the concatenated genome axis. Chromosomes are drawn from position 0 so that
# positions stay comparable to the genome, separated by a small gap.
.genome_layout <- function(d, gap = 0.006) {
  chrs <- intersect(.chr_levels, unique(d$chr))
  # Split the positions once (the slow part on millions of rows), then take
  # each chromosome's range.
  by_chr <- split(d$pos, factor(d$chr, chrs))
  lo <- vapply(by_chr, min, 0, USE.NAMES = FALSE)
  hi <- vapply(by_chr, max, 0, USE.NAMES = FALSE)
  g <- gap * sum(hi)
  offset <- c(0, cumsum(hi + g))[seq_along(chrs)]
  data.frame(chr = chrs, offset = offset, start = offset - g / 2,
             end = offset + hi + g / 2, mid = offset + (lo + hi) / 2,
             index = seq_along(chrs), stringsAsFactors = FALSE)
}

.genome_x <- function(chr, pos, layout) {
  pos + layout$offset[match(chr, layout$chr)]
}

# Drop points that would be drawn on top of an identical-looking point. Only
# points below `below` (in -log10 P) are thinned: the noise floor of a GWAS is
# millions of points that render as a solid band, and keeping one per cell of a
# fine grid changes nothing visible while making plotting fast.
.thin <- function(x, logp, below, nx = 3000, ny_per_unit = 40, group = NULL) {
  low <- logp < below
  rng <- range(x)
  xb <- floor((x - rng[1]) / max(diff(rng), 1) * nx)
  yb <- floor(logp * ny_per_unit)
  key <- xb * 1e6 + yb
  if (!is.null(group)) key <- key * 100 + as.integer(group)
  !(low & duplicated(key))
}

# Compress the y axis above `at` (in -log10 P). Below `at` the axis is
# linear, which is where the decisions are made (the threshold region); above
# it, -log10 P goes on a log scale that gets a fixed share of the panel height
# per decade (`per_decade`, between 12% and 35% in total, and never more than
# a linear axis would give it), so the long tails keep their order without
# eating the plot. Nothing happens unless the data go above the break, and an
# automatic break is only made once the data clear about twice `at`: a break
# for a locus barely past the line (a regional plot whose lead is at 23) is
# clutter for almost no gain.
.y_compressor <- function(max_logp, sig, y_compress = "auto", per_decade = 0.15) {
  at <- if (is.null(y_compress) || isFALSE(y_compress)) {
    NULL
  } else if (is.numeric(y_compress)) {
    y_compress
  } else {
    .auto_compress_at(sig)
  }
  if (is.null(at) || max_logp <= at) {
    return(list(f = identity, at = NULL))
  }
  if (!is.numeric(y_compress) && max_logp <= 2 * at) {
    return(list(f = identity, at = NULL))
  }
  # Never give the tail more height than a linear axis would have.
  frac <- min(0.35, max(0.12, per_decade * log10(max_logp / at)),
              (max_logp - at) / max_logp)
  top <- at / (1 - frac)
  k <- (top - at) / log10(max_logp / at)
  list(f = function(y) ifelse(y <= at, y, at + k * log10(pmax(y, at) / at)),
       at = at)
}

# By default the break sits at about 2.5 times the significance line (20 for
# P = 5e-8), keeping the threshold region large.
.auto_compress_at <- function(sig) {
  .nice_ceiling(max(15, 2.5 * sig))
}

.nice_ceiling <- function(x) {
  step <- if (x <= 20) 5 else if (x <= 100) 10 else if (x <= 500) 50 else 100
  ceiling(x / step) * step
}

.y_breaks <- function(comp, max_logp) {
  if (is.null(comp$at)) {
    b <- scales::breaks_extended(n = 5)(c(0, max_logp))
    b <- b[b >= 0 & b <= max_logp * 1.1]
    return(list(at = b, labels = b))
  }
  lin <- scales::breaks_extended(n = 5)(c(0, comp$at))
  lin <- lin[lin >= 0 & lin <= comp$at]
  # Above the break: 1-2-5 values at least a quarter of a decade apart (on the
  # log scale) from the break and from each other.
  cand <- c(1, 2, 5) * rep(10^(1:4), each = 3)
  cand <- cand[cand > comp$at & cand <= max_logp * 1.05]
  up <- numeric(0)
  last <- comp$at
  for (v in cand) {
    if (log10(v / last) >= 0.25) {
      up <- c(up, v)
      last <- v
    }
  }
  # A tail too short for those still gets one label, so the jump in the axis
  # numbers shows where the scale changes: the largest finer round value in
  # it, at least a tenth of a decade above the break.
  if (!length(up)) {
    fine <- c(1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8) * rep(10^(1:4), each = 10)
    fine <- fine[fine > comp$at & fine <= max_logp & log10(fine / comp$at) >= 0.1]
    if (length(fine)) up <- max(fine)
  }
  b <- c(lin, up)
  list(at = comp$f(b), labels = b)
}

.signed_sqrt <- function() {
  scales::new_transform(
    "signed_sqrt",
    transform = function(x) sign(x) * sqrt(abs(x)),
    inverse = function(x) sign(x) * x^2
  )
}

# Symmetric breaks for an effect axis. On the square-root scale: the largest
# round value inside the limits, about a quarter of it (half-way up on that
# scale) and zero. On the linear scale: ordinary pretty breaks.
.effect_breaks <- function(lim, scale = c("sqrt", "linear"), log_or = FALSE) {
  scale <- match.arg(scale)
  m <- max(abs(lim))
  inside <- function(b) sort(unique(b[b >= -m - 1e-12 & b <= m + 1e-12]))
  if (scale == "linear" && !log_or) {
    return(inside(scales::breaks_extended(n = 5)(c(-m, m))))
  }
  if (log_or) {
    tops <- log(c(1.02, 1.05, 1.1, 1.2, 1.5, 2, 3, 5, 10, 20, 50))
    mids <- tops
  } else {
    tops <- c(1, 2, 3, 4, 5) * rep(10^(-4:3), each = 5)
    mids <- c(1, 2, 5) * rep(10^(-4:3), each = 3)
  }
  top <- max(c(tops[tops <= m], tops[1]))
  mid <- mids[which.min(abs(log(mids / (top / if (scale == "sqrt") 4 else 2))))]
  if (mid >= top) mid <- NULL
  inside(c(-top, -mid, 0, mid, top))
}

# Tick labels for the effect axis, as plain ASCII text ("-0.1"). They are drawn
# through plotmath (see .as_plotmath), which sets a true minus sign on every
# graphics device, including pdf().
.effect_labels <- function(b, log_or = FALSE) {
  v <- if (log_or) exp(b) else b
  out <- trimws(formatC(v, format = "fg", digits = 2))
  out[abs(b) < 1e-12] <- if (log_or) "1" else "0"
  out
}

.as_plotmath <- function(x) {
  x[is.na(x)] <- "''"
  parse(text = x)
}

# The effect axis. "auto": linear, so effects are drawn in proportion, unless
# a few estimates dwarf the rest (more than twice a break set at 3x the median
# lead effect); then the axis breaks there and the outliers go on a log scale
# in the top 25-45% of the track, as on the P-value axis. "sqrt": a smooth
# square-root scale. "linear": no break.
.effect_axis <- function(lead_abs, top, scale = c("auto", "sqrt", "linear"),
                         max_estimate = max(lead_abs)) {
  scale <- match.arg(scale)
  linear <- list(trans = scales::transform_identity(), at = NULL, scale = "linear")
  if (scale == "sqrt") return(list(trans = .signed_sqrt(), at = NULL, scale = "sqrt"))
  if (scale == "linear" || !length(lead_abs)) return(linear)
  at <- .nice_up(3 * stats::median(lead_abs))
  # Break for outlying estimates only; wide intervals alone don't count.
  if (!is.finite(at) || max_estimate <= 2 * at) return(linear)
  decades <- log10(top / at)
  frac <- min(0.45, max(0.25, 0.35 * decades), (top - at) / top)
  k <- (at / (1 - frac) - at) / decades
  list(trans = .broken_log_transform(at, k), at = at, scale = "broken")
}

.nice_up <- function(x) {
  nice <- c(1, 2, 3, 5) * rep(10^(-5:4), each = 4)
  min(nice[nice >= x])
}

# Linear up to `at` (in absolute value), log10 above it; symmetric about zero.
.broken_log_transform <- function(at, k) {
  scales::new_transform(
    "broken_log",
    transform = function(x) {
      a <- abs(x)
      sign(x) * ifelse(a <= at, a, at + k * log10(pmax(a, at) / at))
    },
    inverse = function(u) {
      a <- abs(u)
      sign(u) * ifelse(a <= at, a, at * 10^((a - at) / k))
    }
  )
}

# Ticks for the effect axis.
.track_breaks <- function(ax, top, signed, log_or, roomy = FALSE) {
  if (is.null(ax$at)) {
    if (ax$scale == "linear" && !signed && !log_or) {
      b <- scales::breaks_extended(n = 4)(c(0, top))
      return(b[b >= 0 & b <= top + 1e-12])
    }
    return(.effect_breaks(c(-top, top), if (ax$scale == "sqrt") "sqrt" else "linear", log_or))
  }
  if (log_or) {
    cand <- log(c(1.02, 1.05, 1.1, 1.2, 1.5, 2, 3, 5, 10, 20, 50, 100))
    lin <- utils::tail(cand[cand <= ax$at], 2)
  } else {
    cand <- c(1, 2, 5) * rep(10^(-4:3), each = 3)
    lin <- scales::breaks_extended(n = 3)(c(0, ax$at))
    lin <- lin[lin > 0 & lin <= ax$at]
  }
  up <- numeric(0)
  last <- ax$at
  # The compressed zone is only part of a short track: keep its ticks at least
  # ~half a decade apart.
  for (v in cand[cand > ax$at & cand <= top * 1.05]) {
    if (log10(v / last) >= 0.45) {
      up <- c(up, v)
      last <- v
    }
  }
  # A compressed zone without a tick can't be read at all: fall back to the
  # largest round value inside it, clear of the break.
  if (!length(up)) {
    fine <- if (log_or) {
      log(c(1.1, 1.2, 1.3, 1.5, 1.75, 2, 2.5, 3, 4, 5, 7.5, 10, 15, 20, 30, 50, 75, 100))
    } else {
      c(1, 1.5, 2, 2.5, 3, 4, 5, 7) * rep(10^(-4:3), each = 8)
    }
    tr <- ax$trans$transform
    ok <- fine[fine > ax$at & fine <= top & tr(fine) - tr(ax$at) >= 0.1 * tr(top)]
    up <- utils::tail(ok, 1)
  }
  # On a signed track each half is short: drop the tick at the break itself
  # when there is one above it, or the two collide. A full-height (or
  # full-width) panel has room for both, and the break is the tick that says
  # where the scale changes.
  if (signed && !roomy && length(up)) lin <- lin[abs(lin - ax$at) > 1e-12]
  b <- c(0, lin, up)
  if (signed) b <- c(-rev(b[b > 0]), b)
  b
}
