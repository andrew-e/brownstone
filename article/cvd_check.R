# Colour-vision-deficiency check of the colour sets quoted in article/paper.md:
# smallest pairwise colour distance (delta-E, OKLab x 100) after simulating
# protanopia, deuteranopia and tritanopia (Machado et al. 2009, via colorspace).
# LocusZoom.js values: its default LD bins (statgen/locuszoom,
# esm/layouts/index.js, master branch, read 2026-10-28).
# Usage: Rscript article/cvd_check.R

suppressPackageStartupMessages(library(colorspace))
rgb2hex <- function(m) rgb(m[1], m[2], m[3], maxColorValue = 255)
lz <- c(rgb2hex(c(70,54,153)), rgb2hex(c(38,188,225)), rgb2hex(c(110,254,104)), rgb2hex(c(248,195,42)), rgb2hex(c(219,61,17)))
lz_all <- c(lz, "#9632b8", "#AAAAAA")                       # + lead, + no LD data
names(lz_all) <- c("r2<0.2","0.2-0.4","0.4-0.6","0.6-0.8",">=0.8","lead","no LD")
ours <- c("#c1a8f0", "#9974d4", "#7044ac", "#471d77", "#1e043a")
names(ours) <- c("<0.2","0.2-0.4","0.4-0.6","0.6-0.8",">=0.8")
maf <- c("#2a78d6","#1baf7a","#eb6834","#8f2a5e"); names(maf) <- c(">=5%","1-5%","0.1-1%","<0.1%")
sim <- list(normal = identity, protan = protan, deutan = deutan, tritan = tritan)
oklab <- function(hex) farver::convert_colour(t(col2rgb(hex)), "rgb", "oklab")
dE <- function(hex) { o <- oklab(hex); as.matrix(dist(o)) * 100 }
report <- function(label, pal, ordered = FALSE) {
  cat("\n==", label, "\n")
  for (nm in names(sim)) {
    h <- sim[[nm]](pal); d <- dE(h); diag(d) <- NA
    i <- which(d == min(d, na.rm = TRUE), arr.ind = TRUE)[1, ]
    adj <- if (ordered) min(d[cbind(1:(length(pal)-1), 2:length(pal))]) else NA
    L <- oklab(h)[, 1]
    cat(sprintf("%-7s min pairwise dE %5.1f (%s vs %s)%s  L monotone: %s\n", nm, min(d, na.rm = TRUE),
        names(pal)[i[1]], names(pal)[i[2]],
        if (ordered) sprintf("; min adjacent %5.1f", adj) else "",
        if (ordered) (all(diff(L) > 0) || all(diff(L) < 0)) else "n/a"))
  }
}
report("LocusZoom LD bins (5, ordered)", lz_all[1:5], TRUE)
report("LocusZoom bins + lead + no-data", lz_all)
report("maffect r2 ramp (5, ordered)", ours, TRUE)
report("maffect MAF classes", maf)

# locuszoomr 1.2.0 default LD_scheme (myles-lewis/locuszoomr, R/scatter_plot.R, main):
# grey (no LD), royalblue, cyan2, green3, orange, red (r2 bins), purple (lead).
lzr <- setNames(rgb(t(col2rgb(c("royalblue", "cyan2", "green3", "orange", "red"))), maxColorValue = 255),
                c("<0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", ">=0.8"))
report("locuszoomr LD bins (5, ordered)", lzr, TRUE)
lzr_all <- c(lzr, lead = "#A020F0", `no LD` = "#BEBEBE")
report("locuszoomr bins + lead + no-data", lzr_all)

# Our r2 ramp against our frequency colours (they appear in different plots, but
# a reader moves between them): smallest distance between any ramp bin and any class.
cross <- function(a, b) {
  for (nm in names(sim)) {
    ha <- sim[[nm]](a); hb <- sim[[nm]](b)
    d <- as.matrix(dist(rbind(oklab(ha), oklab(hb))))[seq_along(a), length(a) + seq_along(b)] * 100
    i <- which(d == min(d), arr.ind = TRUE)[1, ]
    cat(sprintf("%-7s min ramp-vs-class dE %5.1f (r2 %s vs MAF %s)\n", nm, min(d), names(a)[i[1]], names(b)[i[2]]))
  }
}
cat("\n== maffect r2 ramp vs maffect MAF classes\n"); cross(ours, maf)
