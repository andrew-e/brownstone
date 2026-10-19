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

# Genes overlapping chr:start-end, from the bundled table for `build` or a
# custom table (columns chr, start, end, gene and optionally strand).
.genes_in_window <- function(chr, start, end, build, genes = NULL) {
  g <- genes %||% .gene_table(build)
  g <- g[.norm_chr(g$chr) == chr & g$end >= start & g$start <= end, , drop = FALSE]
  data.frame(gene = as.character(g$gene), start = as.numeric(g$start), end = as.numeric(g$end),
             strand = if (is.null(g$strand)) rep("+", nrow(g)) else as.character(g$strand),
             stringsAsFactors = FALSE)
}

# Rows for a gene track: first fit, the gene at the lead variant first, then
# by distance from it, so the genes nearest the signal always get a place.
# Each gene takes the span of its body and of its name below it; the name's
# width is estimated (in base pairs) for a panel `width_in` inches wide, as
# text can't be measured before drawing. Genes beyond `max_rows` rows get
# row NA. Returned in that order (nearest first), with each name's position
# and alignment.
.pack_genes <- function(g, start, end, lead_pos, base_size = 10, max_rows = 6, width_in = 7) {
  if (!nrow(g)) {
    g$hjust <- g$label_x <- numeric(0)
    g$row <- integer(0)
    return(g)
  }
  bp_per_in <- (end - start) / width_in
  half <- nchar(g$gene) * 0.55 * (.text_size$note * base_size / 72) * bp_per_in * 1.2 / 2
  mid <- (pmax(g$start, start) + pmin(g$end, end)) / 2
  # A name that would cross an edge is aligned to it; the rest are centred on
  # their gene.
  g$hjust <- ifelse(mid < start + half, 0, ifelse(mid > end - half, 1, 0.5))
  g$label_x <- ifelse(g$hjust == 0, start, ifelse(g$hjust == 1, end, mid))
  lo <- pmin(g$start, g$label_x - 2 * half * g$hjust)
  hi <- pmax(g$end, g$label_x + 2 * half * (1 - g$hjust))
  pad <- 0.02 * (end - start)
  dist <- pmax(g$start - lead_pos, lead_pos - g$end, 0)
  g$row <- NA_integer_
  rows <- list()   # the occupied spans of each row
  ord <- order(dist, g$start)
  for (k in ord) {
    fits <- vapply(rows, function(r) all(lo[k] > r[, 2] + pad | hi[k] < r[, 1] - pad), logical(1))
    r <- which(fits)[1]
    if (is.na(r)) {
      if (length(rows) >= max_rows) next
      rows[[length(rows) + 1L]] <- matrix(c(lo[k], hi[k]), ncol = 2)
      r <- length(rows)
    } else {
      rows[[r]] <- rbind(rows[[r]], c(lo[k], hi[k]))
    }
    g$row[k] <- r
  }
  g[ord, ]
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
