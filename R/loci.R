#' Find independent loci by distance-based clumping
#'
#' Takes the most significant variant, claims every significant variant within
#' `window` base pairs of it as one locus, and repeats until no significant
#' variants are left. This is the usual quick way of turning a GWAS into a list
#' of "hits" for labelling; it does not use LD, so two independent signals
#' closer than `window` are reported as one locus.
#'
#' @param data GWAS summary statistics (anything [standardise_sumstats()] accepts).
#' @param threshold P-value a lead variant must reach.
#' @param window Half-width, in base pairs, of the region each lead variant
#'   claims.
#' @param annotate Gene annotation for each lead variant: `"nearest_gene"` (the
#'   default) uses the bundled Ensembl protein-coding genes for `build`; `NULL`
#'   or `FALSE` skips annotation.
#' @param build Genome build of the positions, `"GRCh38"` or `"GRCh37"`.
#' @param ... Passed on to [standardise_sumstats()] (e.g. column names).
#'
#' @return A data frame with one row per locus, most significant first: the
#'   lead variant (`chr`, `pos`, `snp`, `logp`, `p`, `beta`, `se`, `maf`), the
#'   span of significant variants it claimed (`start`, `end`, `n_sig`), and,
#'   when annotated, `gene` and `gene_distance` (0 when the lead variant is
#'   inside the gene).
#' @export
#' @examples
#' gwas <- simulate_gwas()
#' find_loci(gwas)
find_loci <- function(data, threshold = 5e-8, window = 5e5,
                      annotate = "nearest_gene",
                      build = c("GRCh38", "GRCh37"), ...) {
  build <- match.arg(build)
  d <- .standardise(data, ..., ids = FALSE)
  cl <- .clump(d, -log10(threshold), window)
  cols <- c("chr", "pos", "snp", "logp", "beta", "se", "maf")
  n <- length(cl$lead)
  out <- data.frame(locus = seq_len(n), d[cl$lead, cols, drop = FALSE],
                    stringsAsFactors = FALSE)
  out$snp <- .variant_id(d, cl$lead)
  out$p <- 10^(-out$logp)
  member_pos <- d$pos[cl$member]
  out$start <- as.numeric(vapply(seq_len(n), function(k) min(member_pos[cl$locus_of == k]), 0))
  out$end <- as.numeric(vapply(seq_len(n), function(k) max(member_pos[cl$locus_of == k]), 0))
  out$n_sig <- tabulate(cl$locus_of, nbins = n)
  rownames(out) <- NULL
  if (!n) return(out)

  if (isTRUE(annotate == "nearest_gene")) {
    ng <- nearest_gene(out$chr, out$pos, build = build)
    out$gene <- ng$gene
    out$gene_distance <- ng$distance
  }
  out
}
