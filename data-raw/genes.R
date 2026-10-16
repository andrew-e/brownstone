# Build the bundled gene coordinate tables (genes_grch37, genes_grch38).
#
# Source: Ensembl GTF annotations
#   GRCh37: https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz
#   GRCh38: https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.chr.gtf.gz
#
# Usage (from the package root):
#   Rscript data-raw/genes.R <path/to/GRCh37.gtf.gz> <path/to/GRCh38.gtf.gz>
#
# Only named protein-coding genes on the primary assembly (1-22, X, Y, MT) are
# kept: these are what people almost always want to see as a Manhattan label.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2, all(file.exists(args)))

read_gtf_genes <- function(path, release) {
  cmd <- sprintf("gunzip -c %s | awk -F'\\t' '$3 == \"gene\"'", shQuote(path))
  lines <- readLines(pipe(cmd))
  f <- strsplit(lines, "\t", fixed = TRUE)
  attr_of <- function(x, key) {
    m <- regmatches(x, regexpr(paste0(key, ' "[^"]*"'), x))
    out <- rep(NA_character_, length(x))
    hit <- grepl(paste0(key, ' "'), x)
    out[hit] <- sub(paste0(key, ' "([^"]*)"'), "\\1", m)
    out
  }
  attrs <- vapply(f, `[`, "", 9)
  genes <- data.frame(
    chr = vapply(f, `[`, "", 1),
    start = as.integer(vapply(f, `[`, "", 4)),
    end = as.integer(vapply(f, `[`, "", 5)),
    strand = vapply(f, `[`, "", 7),
    gene = attr_of(attrs, "gene_name"),
    gene_id = attr_of(attrs, "gene_id"),
    biotype = attr_of(attrs, "gene_biotype"),
    stringsAsFactors = FALSE
  )
  # Clone-based placeholder names (RP11-10A14.4, AL627309.1, ...) are never
  # what anyone wants as a label; GRCh37 has ~900 of them.
  clone <- grepl(paste0("^(RP[0-9]+|CT[A-Z]|CTC|CTD|KB|XX|AC[0-9]+|AL[0-9]+|",
                        "AP[0-9]+|AF[0-9]+|BX[0-9]+|CR[0-9]+|Z[0-9]+|U[0-9]+)[-.]|",
                        "^A[A-Z][0-9]{6}\\.[0-9]+$"), genes$gene)
  keep <- genes$biotype == "protein_coding" &
    !is.na(genes$gene) & genes$gene != "" & !clone &
    genes$chr %in% c(1:22, "X", "Y", "MT")
  genes <- genes[keep, c("chr", "start", "end", "strand", "gene", "gene_id")]
  chr_order <- c(1:22, "X", "Y", "MT")
  genes <- genes[order(match(genes$chr, chr_order), genes$start), ]
  rownames(genes) <- NULL
  attr(genes, "source") <- release
  genes
}

genes_grch37 <- read_gtf_genes(args[1], "Ensembl GRCh37 release 87")
genes_grch38 <- read_gtf_genes(args[2], "Ensembl release 116 (GRCh38)")

message("GRCh37: ", nrow(genes_grch37), " genes; GRCh38: ", nrow(genes_grch38), " genes")

dir.create("data", showWarnings = FALSE)
save(genes_grch37, file = "data/genes_grch37.rda", compress = "xz")
save(genes_grch38, file = "data/genes_grch38.rda", compress = "xz")
