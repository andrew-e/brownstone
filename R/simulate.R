.chr_lengths <- list(
  GRCh38 = c(`1` = 248956422, `2` = 242193529, `3` = 198295559, `4` = 190214555,
             `5` = 181538259, `6` = 170805979, `7` = 159345973, `8` = 145138636,
             `9` = 138394717, `10` = 133797422, `11` = 135086622, `12` = 133275309,
             `13` = 114364328, `14` = 107043718, `15` = 101991189, `16` = 90338345,
             `17` = 83257441, `18` = 80373285, `19` = 58617616, `20` = 64444167,
             `21` = 46709983, `22` = 50818468, X = 156040895, Y = 57227415),
  GRCh37 = c(`1` = 249250621, `2` = 243199373, `3` = 198022430, `4` = 191154276,
             `5` = 180915260, `6` = 171115067, `7` = 159138663, `8` = 146364022,
             `9` = 141213431, `10` = 135534747, `11` = 135006516, `12` = 133851895,
             `13` = 115169878, `14` = 107349540, `15` = 102531392, `16` = 90354753,
             `17` = 81195210, `18` = 78077248, `19` = 59128983, `20` = 63025520,
             `21` = 48129895, `22` = 51304566, X = 155270560, Y = 59373566)
)

#' Causal loci used by [simulate_gwas()]
#'
#' A mix of common, low-frequency and rare causal variants, each placed inside
#' a real gene so that gene annotation can be demonstrated. The genes are ones
#' with well-known lipid associations, but **the effect sizes and P-values are
#' invented**: they only mimic the shape of a real GWAS (common variants with
#' small effects, rare variants with large ones).
#'
#' @format A data frame with columns `gene`, `maf` (minor allele frequency of
#'   the causal variant) and `z` (its z statistic for the minor allele).
#' @export
sim_loci <- data.frame(
  gene = c(
    # common
    "CELSR2", "APOB", "ABCG8", "HMGCR", "TRIB1", "LDLR", "APOE", "HNF1A",
    "ABO", "FADS2", "MAFB", "TOP1", "LPA", "GCKR", "NPC1L1", "PPP1R3B",
    "TIMD4", "LIPC", "CETP", "ANGPTL3", "ST3GAL4", "MYLIP", "FRK",
    # low frequency
    "PCSK9", "TM6SF2", "LPL",
    # rare
    "LDLR", "APOC3", "ABCA1", "LIPG", "ASGR1", "SCARB1"),
  maf = c(0.22, 0.30, 0.07, 0.38, 0.45, 0.11, 0.08, 0.35, 0.20, 0.33, 0.30,
          0.15, 0.07, 0.40, 0.25, 0.15, 0.30, 0.35, 0.32, 0.32, 0.25, 0.28,
          0.20,
          0.015, 0.035, 0.03,
          0.002, 0.0007, 0.004, 0.006, 0.0008, 0.005),
  z = c(-18, 14, -10, 11, -8.5, -13, -38, 7, 9, 8, 6.2, 5.8, 12, 6.5, 5.9,
        6.8, 7.5, 6, 9.5, 7.2, 5.2, 5.4, 4.9,
        -14, -9, 7,
        12, -9, -7.5, 6.8, -8.2, 6.5),
  stringsAsFactors = FALSE
)

#' Simulate GWAS summary statistics
#'
#' Simulates a genome-wide association study of a quantitative trait with
#' realistic-looking structure: an allele frequency spectrum rich in rare
#' variants, standard errors that grow as variants get rarer, and LD around each
#' causal variant. Because a rare allele cannot be in strong LD with a common
#' one, rare causal variants stand alone while common ones raise a tower of
#' correlated neighbours, as in real data.
#'
#' @param n_variants Number of variants.
#' @param n Sample size (sets the standard errors).
#' @param build Genome build for positions, `"GRCh38"` or `"GRCh37"`.
#' @param loci Causal loci; see [sim_loci].
#' @param n_background Number of weak (sub-threshold) signals to scatter
#'   across the genome.
#' @param inflation Variance of the null z statistics (a little above 1, as in a
#'   real polygenic trait).
#' @param seed Random seed, or `NULL` to use the current RNG state. The global
#'   RNG state is restored afterwards.
#'
#' @return A data frame with the columns of a typical summary statistics file:
#'   `CHR`, `BP`, `SNP`, `A1` (effect allele), `A2`, `EAF` (effect allele
#'   frequency), `BETA`, `SE` and `P`. The causal variants are attached as the
#'   attribute `"causal"`.
#' @export
#' @examples
#' gwas <- simulate_gwas(n_variants = 1e5)
#' head(gwas)
simulate_gwas <- function(n_variants = 1e6, n = 150000,
                          build = c("GRCh38", "GRCh37"), loci = sim_loci,
                          n_background = 40, inflation = 1.04, seed = 42) {
  build <- match.arg(build)
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
      on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
  }
  chr_len <- .chr_lengths[[build]][as.character(1:22)]
  n_chr <- as.vector(stats::rmultinom(1, n_variants, chr_len))
  chr <- rep(1:22, n_chr)
  pos <- unlist(lapply(1:22, function(i) {
    sort(sample.int(chr_len[[i]] - 200000L, n_chr[i]) + 100000L)
  }))

  # Allele frequency spectrum of a well-imputed or sequenced GWAS, down to
  # MAF ~ 0.025% (ultra-rare < 0.1%, rare 0.1-1%, low-frequency 1-5%, common).
  u <- stats::runif(n_variants)
  maf <- ifelse(u < 0.12, 10^stats::runif(n_variants, -3.6, -3),
         ifelse(u < 0.40, 10^stats::runif(n_variants, -3, -2),
         ifelse(u < 0.55, 10^stats::runif(n_variants, -2, log10(0.05)),
                stats::runif(n_variants, 0.05, 0.5))))
  z <- stats::rnorm(n_variants, sd = sqrt(inflation))

  add_signal <- function(c_idx, zc, ld_scale) {
    on <- which(chr == chr[c_idx])
    win <- on[abs(pos[on] - pos[c_idx]) <= 1e6]
    win <- win[win != c_idx]
    d <- abs(pos[win] - pos[c_idx])
    # A haplotype block of tight proxies with frequencies close to the causal.
    block <- d < 1.5 * ld_scale & maf[win] >= 0.05 & stats::runif(length(win)) < 0.6
    if (maf[c_idx] >= 0.05) {
      maf[win[block]] <<- pmin(0.5, maf[c_idx] * exp(stats::rnorm(sum(block), 0, 0.08)))
    }
    f1 <- pmin(maf[win], maf[c_idx])
    f2 <- pmax(maf[win], maf[c_idx])
    r_max <- sqrt(f1 * (1 - f2) / (f2 * (1 - f1)))
    r <- exp(-d / ld_scale) * stats::rbeta(length(win), 1.6, 1)
    r[block] <- stats::runif(sum(block), 0.6, 1)
    r <- r * r_max * ifelse(stats::runif(length(win)) < 0.15, -1, 1)
    z[win] <<- z[win] * sqrt(1 - r^2) + r * zc
    z[c_idx] <<- zc + stats::rnorm(1, 0, 0.3)
  }

  # Causal loci, placed inside their genes.
  genes <- .gene_table(build)
  causal <- data.frame(gene = loci$gene, chr = rep(NA_character_, nrow(loci)),
                       pos = rep(NA_real_, nrow(loci)), maf = loci$maf, z = loci$z,
                       stringsAsFactors = FALSE)
  for (k in seq_len(nrow(loci))) {
    g <- genes[genes$gene == loci$gene[k], , drop = FALSE][1, ]
    if (is.na(g$chr) || !g$chr %in% as.character(1:22)) next
    on <- which(chr == as.integer(g$chr))
    target <- g$start + stats::runif(1) * (g$end - g$start)
    c_idx <- on[which.min(abs(pos[on] - target))]
    maf[c_idx] <- loci$maf[k]
    add_signal(c_idx, loci$z[k], stats::runif(1, 20e3, 100e3))
    causal$chr[k] <- g$chr
    causal$pos[k] <- pos[c_idx]
  }

  # Weak background signals, plus two lone rare variants of the kind that
  # usually turn out to be imputation artefacts.
  for (k in seq_len(n_background)) {
    c_idx <- sample.int(n_variants, 1)
    maf[c_idx] <- stats::runif(1, 0.05, 0.5)
    add_signal(c_idx, sample(c(-1, 1), 1) * stats::runif(1, 3.2, 5.0),
               stats::runif(1, 15e3, 60e3))
  }
  for (k in 1:2) {
    c_idx <- sample.int(n_variants, 1)
    maf[c_idx] <- c(0.0004, 0.0015)[k]
    z[c_idx] <- sample(c(-1, 1), 1) * stats::runif(1, 5.6, 6.2)
  }

  se <- 1 / sqrt(2 * n * maf * (1 - maf)) * exp(stats::rnorm(n_variants, 0, 0.03))
  # Report most variants with the minor allele as the effect allele, as most
  # pipelines do (effect allele = ALT, usually the minor allele).
  minor_is_effect <- stats::runif(n_variants) < 0.8
  z_eff <- ifelse(minor_is_effect, z, -z)
  bases <- c("A", "C", "G", "T")
  a1 <- sample(bases, n_variants, replace = TRUE)
  a2 <- bases[(match(a1, bases) + sample.int(3, n_variants, replace = TRUE) - 1L) %% 4L + 1L]

  out <- data.frame(
    CHR = chr,
    BP = pos,
    SNP = paste(chr, pos, a2, a1, sep = ":"),
    A1 = a1,
    A2 = a2,
    EAF = signif(ifelse(minor_is_effect, maf, 1 - maf), 4),
    BETA = signif(z_eff * se, 5),
    SE = signif(se, 5),
    P = signif(2 * stats::pnorm(-abs(z_eff)), 4),
    stringsAsFactors = FALSE
  )
  attr(out, "causal") <- causal
  out
}
