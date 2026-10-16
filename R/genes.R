#' Protein-coding gene coordinates
#'
#' Named protein-coding genes from Ensembl, used to label loci with their
#' nearest gene. Built by `data-raw/genes.R`.
#'
#' @format A data frame with columns `chr`, `start`, `end` (1-based, inclusive),
#'   `strand`, `gene` (HGNC symbol) and `gene_id` (Ensembl ID).
#' @source `genes_grch37`: Ensembl GRCh37 release 87.
#'   `genes_grch38`: Ensembl release 116.
#'   <https://ftp.ensembl.org/pub/>
#' @name genes
"genes_grch37"

#' @rdname genes
"genes_grch38"

.gene_table <- function(build) {
  switch(build,
         GRCh38 = brownstone::genes_grch38,
         GRCh37 = brownstone::genes_grch37,
         stop("Unknown build '", build, "'.", call. = FALSE))
}

#' Nearest gene to each position
#'
#' For each position, returns the gene whose body contains it, or failing that
#' the gene with the closest end. When a position lies inside several genes, the
#' one whose transcription start site is closest wins.
#'
#' @param chr,pos Chromosomes and base-pair positions (recycled together).
#' @param build Genome build of the positions, `"GRCh38"` or `"GRCh37"`.
#' @param genes Optional custom gene table with columns `chr`, `start`, `end`,
#'   `gene` and (optionally) `strand`, overriding `build`.
#' @param max_distance Positions further than this from every gene get `NA`.
#'
#' @return A data frame with columns `gene` and `distance` (base pairs; 0 inside
#'   a gene).
#' @export
#' @examples
#' nearest_gene(c("1", "19"), c(55050000, 44905000))
nearest_gene <- function(chr, pos, build = c("GRCh38", "GRCh37"), genes = NULL,
                         max_distance = Inf) {
  build <- match.arg(build)
  if (is.null(genes)) genes <- .gene_table(build)
  chr <- .norm_chr(chr)
  g_chr <- .norm_chr(genes$chr)
  strand <- if (is.null(genes$strand)) rep("+", nrow(genes)) else genes$strand
  tss <- ifelse(strand == "-", genes$end, genes$start)
  n <- max(length(chr), length(pos))
  chr <- rep_len(chr, n)
  pos <- rep_len(pos, n)
  gene <- rep(NA_character_, n)
  dist <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    on <- which(g_chr == chr[i])
    if (!length(on) || !is.finite(pos[i])) next
    d <- pmax(genes$start[on] - pos[i], pos[i] - genes$end[on], 0)
    best <- on[d == min(d)]
    if (length(best) > 1L) best <- best[which.min(abs(tss[best] - pos[i]))]
    if (min(d) <= max_distance) {
      gene[i] <- genes$gene[best]
      dist[i] <- min(d)
    }
  }
  data.frame(gene = gene, distance = dist, stringsAsFactors = FALSE)
}
