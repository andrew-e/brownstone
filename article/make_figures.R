# Regenerate the figures in article/paper.md.
#
# Usage (from the package root):
#   Rscript article/make_figures.R <path/to/finngen_R12_E4_HYPERCHOL.gz>
#
# The FinnGen release 12 summary statistics for E4_HYPERCHOL are public:
#   https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_E4_HYPERCHOL.gz
# Without the file, figures 1 and 4 (simulated data) are still drawn. If a data
# frame `fg` of the FinnGen statistics already exists, it is used instead of
# reading the file.

suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
out <- "article/figures"
save <- function(p, file, w = 12, h = 6.5) {
  ggplot2::ggsave(file.path(out, file), p, width = w, height = h, dpi = 200,
                  device = ragg::agg_png, bg = "white")
  message("wrote ", file)
}

# Figure 1: simulated GWAS with the suggestive tier.
gwas <- simulate_gwas()
save(manhattan(gwas, suggestive = 1e-5), "figure1_simulated.png")

# Figure 4: the other plots, on simulated data. A: volcano plot. B: regional
# plot of the LDLR locus, with the LD simulated alongside the data. C: Miami
# plot of two separately simulated cohorts; the second lacks some of the first's
# loci, is weaker at others and has loci of its own (as in the README).
women_loci <- sim_loci[!sim_loci$gene %in% c("HNF1A", "TOP1", "MAFB", "LIPG", "ASGR1"), ]
women_loci$z <- women_loci$z * ifelse(women_loci$gene %in% c("APOE", "LPA"), 0.6,
                                       ifelse(women_loci$gene %in% c("ABO", "FADS2", "CELSR2"), 1.3, 1))
women_loci <- rbind(women_loci, data.frame(
  gene = c("PNPLA3", "GPAM", "HFE", "PLTP", "ANGPTL4", "LCAT"),
  maf = c(0.25, 0.30, 0.07, 0.18, 0.004, 0.003),
  z = c(12, 8, 6.5, 6, -9, 8.5)))
women <- simulate_gwas(n = 120000, loci = women_loci, seed = 7)
pad <- ggplot2::theme(plot.margin = ggplot2::margin(t = 22, r = 6, b = 4, l = 22))
suite <- suppressMessages({
  a <- volcano_plot(gwas, max_labels = 12) + pad
  b <- region_plot(gwas, gene = "LDLR", ld = attr(gwas, "ld")) & pad
  c <- miami_plot(gwas, women, names = c("Men", "Women"), max_labels = 12) & pad
  (patchwork::wrap_elements(full = a) | patchwork::wrap_elements(full = b)) /
    patchwork::wrap_elements(full = c) +
    patchwork::plot_layout(heights = c(1, 1.5)) +
    patchwork::plot_annotation(tag_levels = "A")
})
save(suite, "figure4_suite.png", w = 13, h = 16)

# Figures 2 and 3: FinnGen R12 pure hypercholesterolaemia (E4_HYPERCHOL).
args <- commandArgs(trailingOnly = TRUE)
if (!exists("fg") && length(args) && file.exists(args[1])) {
  fg <- data.table::fread(args[1], data.table = FALSE, showProgress = FALSE)
}
if (exists("fg")) {
  save(manhattan(fg, build = "GRCh38"), "figure2_finngen_hyperchol.png")
  save(effect_vs_maf(fg), "figure3_finngen_effect_vs_maf.png", w = 9, h = 6.5)
} else {
  message("FinnGen file not given; skipping figures 2 and 3.")
}
