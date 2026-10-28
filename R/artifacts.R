# The artefact check: variants whose reported effects disagree with their own
# case and control allele frequencies.

#' Variants whose effects disagree with their own allele frequencies
#'
#' @description
#' For a binary trait whose summary statistics report the allele frequency
#' among cases and among controls (FinnGen's `af_alt_cases`/`af_alt_controls`,
#' SAIGE's `AF_case`/`AF_ctrl`), the two frequencies imply a crude odds ratio,
#' and a genuine association should roughly agree with the reported estimate:
#' covariate adjustment moves an odds ratio somewhat, but not by an order of
#' magnitude. `find_artifacts()` screens the lead variant of every
#' genome-wide significant locus and returns the comparison for each.
#'
#' A lead is `flagged` when the reported odds ratio and the crude one differ by
#' at least `ratio`-fold, or the reported effect (on the log scale) differs
#' from the crude one by at least `z` standard errors. The two frequencies are
#' first folded to the same orientation as the reported effect (see
#' `effect_allele` in [standardise_sumstats()]), so the comparison is fair.
#' Frequencies of exactly 0, 1 or missing cannot be judged (0 and 1 may only
#' be rounding) and are left out. `direction_disagrees` is reported as well,
#' but a sign flip alone is not flagged: with a handful of carriers the crude
#' odds ratio is noise.
#'
#' This is a screen, not a test. Its typical catch is an ultra-rare variant
#' whose saddlepoint or Firth test produced an estimate its own allele counts
#' don't support; the numbers are all returned for judgment.
#'
#' @param data GWAS summary statistics (anything [standardise_sumstats()]
#'   accepts).
#' @param threshold P-value a lead variant must reach.
#' @param window Half-width, in base pairs, of the region each lead variant
#'   claims (see [find_loci()]).
#' @param annotate Gene annotation for each lead variant: `"nearest_gene"`
#'   (the default) uses the bundled Ensembl protein-coding genes for `build`;
#'   `NULL` or `FALSE` skips annotation.
#' @param build Genome build of the positions, `"GRCh38"` or `"GRCh37"`.
#' @param ratio Flag when the reported odds ratio and the crude one differ by
#'   at least this factor (either direction).
#' @param z Flag when the reported effect differs from the crude one by at
#'   least this many standard errors.
#' @param ... Passed on to [standardise_sumstats()] (e.g. column names).
#'
#' @return A data frame, most significant lead first, one row per lead variant
#'   with frequencies that can be judged: `snp`, `chr`, `pos`, `gene` (when
#'   annotated), `logp`, `beta`, `se`, `af_case`, `af_ctrl`, `crude_or` (the
#'   odds ratio the two frequencies imply), `discrepancy` (beta minus the crude
#'   log odds ratio), `ratio` (the factor between the two odds ratios, either
#'   direction), `z_discrepancy` (`|discrepancy| / se`), `direction_disagrees`
#'   and `flagged`. Empty, with a message, when the data have no case and
#'   control frequency columns.
#' @export
#' @examples
#' gwas <- simulate_gwas(n_variants = 2e5)
#' head(find_artifacts(gwas))   # no case/control frequencies: nothing to check
find_artifacts <- function(data, threshold = 5e-8, window = 5e5,
                           annotate = "nearest_gene",
                           build = c("GRCh38", "GRCh37"),
                           ratio = 10, z = 6, ...) {
  build <- match.arg(build)
  d <- .standardise(data, ..., ccafs = TRUE)
  if (!isTRUE(attr(d, "has_ccaf"))) {
    message("No case and control allele frequency columns found (FinnGen's ",
            "af_alt_cases/af_alt_controls, SAIGE's AF_case/AF_ctrl); nothing to check.")
    return(invisible(.artifact_frame()))
  }
  out <- .artifact_table(d, threshold, window, ratio, z)
  rows <- out$row
  out$row <- NULL
  if (!nrow(out)) return(out)
  if (isTRUE(annotate == "nearest_gene")) {
    ng <- nearest_gene(out$chr, out$pos, build = build)
    out$gene <- ng$gene
    out <- out[, c("snp", "chr", "pos", "gene", "logp", "beta", "se",
                   "af_case", "af_ctrl", "crude_or", "discrepancy", "ratio",
                   "z_discrepancy", "direction_disagrees", "flagged")]
  }
  if (any(out$flagged)) {
    .artifact_message(out, build)
  }
  out
}

# The empty artefact table, with all the columns of a full one.
.artifact_frame <- function() {
  data.frame(snp = character(0), chr = character(0), pos = numeric(0),
             logp = numeric(0), beta = numeric(0), se = numeric(0),
             af_case = numeric(0), af_ctrl = numeric(0), crude_or = numeric(0),
             discrepancy = numeric(0), ratio = numeric(0),
             z_discrepancy = numeric(0), direction_disagrees = logical(0),
             flagged = logical(0), row = integer(0), stringsAsFactors = FALSE)
}

# The artefact comparison for the lead variant of every genome-wide
# significant locus in `d`, with `row` (the row of `d`) kept so callers can
# mark the variant on a plot.
.artifact_table <- function(d, threshold, window, ratio = 10, z = 6) {
  cl <- .clump(d, -log10(threshold), window)
  lead <- cl$lead
  af1 <- d$af_case[lead]
  af2 <- d$af_ctrl[lead]
  # A frequency of 0 or 1 may only be rounding: it cannot be judged.
  ok <- is.finite(af1) & is.finite(af2) &
    af1 > 0 & af1 < 1 & af2 > 0 & af2 < 1
  out <- .artifact_frame()
  if (!any(ok)) return(out)
  out <- data.frame(snp = .variant_id(d, lead[ok]), chr = d$chr[lead[ok]],
                    pos = d$pos[lead[ok]], logp = d$logp[lead[ok]],
                    beta = d$beta[lead[ok]], se = d$se[lead[ok]],
                    af_case = af1[ok], af_ctrl = af2[ok],
                    stringsAsFactors = FALSE)
  out$crude_or <- (out$af_case / (1 - out$af_case)) /
    (out$af_ctrl / (1 - out$af_ctrl))
  out$discrepancy <- out$beta - log(out$crude_or)
  out$ratio <- exp(abs(out$discrepancy))
  out$z_discrepancy <- abs(out$discrepancy) / out$se
  out$direction_disagrees <- sign(out$beta) != sign(log(out$crude_or))
  # A sign flip alone is not flagged (see above); an order of magnitude or a
  # 6-SE disagreement is.
  out$flagged <- out$ratio >= ratio |
    (!is.na(out$z_discrepancy) & out$z_discrepancy >= z)
  out$row <- lead[ok]
  out
}

# One message about the flagged variants, shared by find_artifacts() and the
# plots: how many, and the worst (the largest factor between the two odds
# ratios), with its nearest gene.
.artifact_message <- function(out, build) {
  i <- which(out$flagged)[which.max(out$ratio[out$flagged])]
  gene <- nearest_gene(out$chr[i], out$pos[i], build = build)$gene
  gene <- if (is.na(gene)) .chr_pos_id(out$chr[i], out$pos[i]) else gene
  one <- sum(out$flagged) == 1
  message(sum(out$flagged), " lead variant", if (!one) "s",
          " disagree", if (one) "s with its" else " with their",
          " own case/control allele frequencies (worst: ",
          out$snp[i], " near ", gene, "); see find_artifacts() for the numbers.")
}

# The rows, in `d`, of the flagged lead variants, for the plots' flag markers.
# Drawn variants only: called on the thinned data, whose leads all survive
# (hits are never thinned). Quiet when there is nothing to flag.
.artifact_flags <- function(d, threshold, window, build, ratio = 10, z = 6) {
  if (!isTRUE(attr(d, "has_ccaf"))) return(integer(0))
  out <- .artifact_table(d, threshold, window, ratio, z)
  if (!any(out$flagged)) return(integer(0))
  .artifact_message(out, build)
  out$row[out$flagged]
}
