# Column detection and standardisation of GWAS summary statistics.
#
# Every user-facing function accepts summary statistics as they come out of the
# common GWAS tools (PLINK, REGENIE, BOLT-LMM, SAIGE, GWAS-SSF, METAL, ...) and
# guesses the columns. Any guess can be overridden by naming the column.

.aliases <- list(
  chr = c("chr", "chrom", "chromosome", "#chrom", "#chr", "chr_name",
          "chromosome_name", "hg19chrc", "chrom_b37", "chrom_b38"),
  pos = c("pos", "bp", "position", "base_pair_location", "genpos", "bp_hg19",
          "bp_hg38", "pos_b37", "pos_b38", "pos_hg19", "pos_hg38", "basepair",
          "physical_position", "start"),
  p = c("p", "pval", "p_value", "pvalue", "p.value", "p-value", "p_bolt_lmm",
        "p_bolt_lmm_inf", "p_wald", "p_lrt", "frequentist_add_pvalue",
        "all_inv_var_meta_p", "pval_nominal", "gc_pvalue"),
  log10p = c("log10p", "log10_p", "mlogp", "neglog10p", "neglog10_p", "neg_log10_p",
             "mlog10p", "minus_log10_p", "neg_log_10_p_value", "-log10p",
             "-log10(p)", "log10_pval", "lp"),
  beta = c("beta", "b", "effect", "beta_hat", "effect_size", "est", "estimate",
           "frequentist_add_beta_1", "all_inv_var_meta_beta", "log_or", "logor",
           "lnor", "beta_1"),
  or = c("or", "odds_ratio", "oddsratio"),
  se = c("se", "standard_error", "sebeta", "se_beta", "stderr", "std_err",
         "std_error", "frequentist_add_se_1", "all_inv_var_meta_sebeta", "se_1",
         "log_or_se", "se_log_or"),
  maf = c("maf", "minor_allele_frequency", "minor_af"),
  eaf = c("eaf", "effect_allele_frequency", "a1freq", "a1_freq", "freq1",
          "freq", "af", "af_alt", "alt_af", "alt_freq", "af_allele2", "a1_af",
          "freq_tested_allele", "freq_a1", "effect_allele_freq", "frq", "af1",
          "all_meta_af", "frequency"),
  snp = c("rsid", "snp", "rs", "rs_id", "rsids", "markername", "snpid",
          "snp_id", "variant_id", "id", "markerid", "marker", "variant",
          "marker_name")
)

.chr_levels <- c(as.character(1:22), "X", "Y", "XY", "MT")

# Return the column of `data` to use for `field`: the user's choice if given,
# otherwise the first alias that matches (case-insensitively), otherwise NULL.
.find_col <- function(data, field, user = NULL) {
  nms <- names(data)
  if (!is.null(user)) {
    if (!is.character(user) || length(user) != 1L || !user %in% nms) {
      stop(sprintf("Column '%s' (given for `%s`) is not in the data.",
                   paste(user, collapse = ", "), field), call. = FALSE)
    }
    return(user)
  }
  lower <- tolower(nms)
  for (a in .aliases[[field]]) {
    hit <- which(lower == a)
    if (length(hit)) return(nms[hit[1L]])
  }
  NULL
}

.norm_chr <- function(x) {
  # Work on the distinct values only: tens of millions of rows, ~25 values.
  f <- factor(x)
  lev <- toupper(sub("^chr", "", trimws(levels(f)), ignore.case = TRUE))
  lev[lev == "23"] <- "X"
  lev[lev == "24"] <- "Y"
  lev[lev == "25"] <- "XY"
  lev[lev %in% c("26", "M")] <- "MT"
  lev[as.integer(f)]
}

.chr_pos_id <- function(chr, pos) {
  paste(chr, format(pos, scientific = FALSE, trim = TRUE), sep = ":")
}

# Variant IDs for rows `i` of standardised data: the data's own IDs, or, when
# there was no ID column, "chr:pos" built for these rows only.
.variant_id <- function(d, i = seq_len(nrow(d))) {
  if ("snp" %in% names(attr(d, "columns"))) return(d$snp[i])
  .chr_pos_id(d$chr[i], d$pos[i])
}

# -log10(P) from a z statistic, computed on the log scale so that it never
# underflows (a z of 40 is P ~ 1e-350, below the smallest double).
.logp_from_z <- function(z) {
  -(stats::pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10)
}

# -log10(P) from P-values written as text in scientific notation, exact even
# below the smallest double ("1.24e-652" -> 651.9). NA where not parseable.
.logp_from_text <- function(x) {
  x <- trimws(as.character(x))
  ok <- grepl("^[0-9]*\\.?[0-9]+[eE][-+]?[0-9]+$", x)
  out <- rep(NA_real_, length(x))
  mant <- as.numeric(sub("[eE].*$", "", x[ok]))
  expo <- as.numeric(sub("^.*[eE]", "", x[ok]))
  out[ok] <- ifelse(mant > 0, -(log10(mant) + expo), NA_real_)
  out
}

#' Standardise GWAS summary statistics
#'
#' Finds the chromosome, position, P-value, effect size, standard error, allele
#' frequency and variant ID columns in a table of summary statistics, and returns
#' them under fixed names. Column names from PLINK, REGENIE, BOLT-LMM, SAIGE,
#' METAL and the GWAS-SSF standard are recognised automatically; any column can
#' be named explicitly instead.
#'
#' @param data A data frame of GWAS summary statistics, one row per variant.
#' @param chr,pos,p,beta,se,snp Optional column names (strings) overriding the
#'   automatic detection. `p` may also be a column of -log10(P) values if it is
#'   named like one (e.g. `LOG10P`); use `log10p` to say so explicitly.
#' @param log10p Optional name of a column holding -log10(P).
#' @param or Optional name of an odds-ratio column, used (on the log scale) when
#'   there is no beta column.
#' @param maf Optional name of a minor allele frequency column.
#' @param eaf Optional name of an effect (tested) allele frequency column.
#' @param effect_allele How to orient effects. `"minor"` (the default) flips the
#'   sign of beta wherever the effect allele is the major allele, so every effect
#'   reads "per copy of the minor allele" and has the same meaning across the
#'   plot. This needs an effect allele frequency column; without one effects are
#'   used as given. `"as_is"` never flips.
#' @param keep Names of any other columns to carry through unchanged (e.g. a
#'   column of gene names to label hits with).
#' @param quiet If `TRUE`, don't report which columns were used.
#'
#' @return A data frame with columns `chr` (character: 1-22, X, Y, XY, MT),
#'   `pos`, `snp`, `logp` (-log10 P), `beta`, `se`, `maf` and any `keep`
#'   columns, plus the attribute `"columns"` recording which input columns were
#'   used. `beta`, `se` and `maf` are `NA` when the input has no such column.
#' @export
#' @examples
#' gwas <- simulate_gwas(n_variants = 2e4)
#' head(standardise_sumstats(gwas))
standardise_sumstats <- function(data, chr = NULL, pos = NULL, p = NULL,
                                 beta = NULL, se = NULL, snp = NULL,
                                 log10p = NULL, or = NULL, maf = NULL,
                                 eaf = NULL,
                                 effect_allele = c("minor", "as_is"),
                                 keep = NULL, quiet = TRUE) {
  .standardise(data, chr = chr, pos = pos, p = p, beta = beta, se = se, snp = snp,
               log10p = log10p, or = or, maf = maf, eaf = eaf,
               effect_allele = effect_allele, keep = keep, quiet = quiet)
}

# The work of standardise_sumstats(). With `ids = FALSE`, variants are not
# named "chr:pos" when the data have no ID column: building those names for
# millions of rows is slow, and only labelled loci need them (.variant_id()).
.standardise <- function(data, chr = NULL, pos = NULL, p = NULL,
                         beta = NULL, se = NULL, snp = NULL,
                         log10p = NULL, or = NULL, maf = NULL,
                         eaf = NULL,
                         effect_allele = c("minor", "as_is"),
                         keep = NULL, quiet = TRUE, ids = TRUE) {
  if (inherits(data, "brownstone_sumstats") && all(keep %in% names(data))) return(data)
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  effect_allele <- match.arg(effect_allele)
  missing_keep <- setdiff(keep, names(data))
  if (length(missing_keep)) {
    stop("Column(s) not in the data: ", paste(missing_keep, collapse = ", "), call. = FALSE)
  }

  col <- list(
    chr = .find_col(data, "chr", chr),
    pos = .find_col(data, "pos", pos),
    p = .find_col(data, "p", p),
    log10p = .find_col(data, "log10p", log10p),
    beta = .find_col(data, "beta", beta),
    or = .find_col(data, "or", or),
    se = .find_col(data, "se", se),
    maf = .find_col(data, "maf", maf),
    eaf = .find_col(data, "eaf", eaf),
    snp = .find_col(data, "snp", snp)
  )
  if (is.null(col$chr) || is.null(col$pos)) {
    stop("Could not find chromosome and position columns; name them with ",
         "`chr =` and `pos =`.", call. = FALSE)
  }
  if (is.null(col$p) && is.null(col$log10p)) {
    stop("Could not find a P-value column; name it with `p =` (or `log10p =`).",
         call. = FALSE)
  }
  # REGENIE writes one row per test (ADD, DOM, REC, ...) when asked for more
  # than the additive model; mixing them would plot each variant several times.
  test_col <- names(data)[tolower(names(data)) == "test"]
  if (length(test_col) == 1L) {
    tests <- unique(data[[test_col]])
    if (length(tests) > 1L && "ADD" %in% tests) {
      if (!quiet) message("Keeping the additive (TEST == \"ADD\") rows of ", length(tests), " tests.")
      data <- data[data[[test_col]] %in% "ADD", , drop = FALSE]
    }
  }
  n <- nrow(data)
  num <- function(name) if (is.null(name)) rep(NA_real_, n) else as.numeric(data[[name]])

  out <- data.frame(
    chr = .norm_chr(data[[col$chr]]),
    pos = num(col$pos),
    stringsAsFactors = FALSE
  )
  out$snp <- if (!is.null(col$snp)) {
    as.character(data[[col$snp]])
  } else if (ids) {
    .chr_pos_id(out$chr, out$pos)
  } else {
    rep(NA_character_, nrow(out))
  }

  # Effect size: beta, or log(OR).
  out$beta <- if (!is.null(col$beta)) num(col$beta) else if (!is.null(col$or)) log(num(col$or)) else NA_real_
  out$se <- num(col$se)

  # Allele frequency: MAF directly, or folded from the effect allele frequency.
  eaf_v <- num(col$eaf)
  out$maf <- if (!is.null(col$maf)) num(col$maf) else eaf_v
  out$maf <- pmin(out$maf, 1 - out$maf)
  if (effect_allele == "minor" && !is.null(col$eaf)) {
    flip <- !is.na(eaf_v) & eaf_v > 0.5
    out$beta[flip] <- -out$beta[flip]
  }

  # -log10(P). Prefer an explicit -log10 column; otherwise transform P. P-values
  # too small for a double (they read as zero) are recovered exactly from their
  # text if the column was read as text, and otherwise from beta / se.
  if (!is.null(col$log10p)) {
    out$logp <- num(col$log10p)
    # A column of log10(P) (mostly negative) rather than -log10(P). Counting
    # is much faster than a median on millions of rows.
    if (mean(out$logp < 0, na.rm = TRUE) > 0.5) {
      out$logp <- -out$logp
      if (!quiet) message("Column '", col$log10p, "' holds log10(P); using its negative.")
    }
  } else {
    raw <- data[[col$p]]
    pv <- if (is.numeric(raw)) raw else suppressWarnings(as.numeric(as.character(raw)))
    if (any(pv < 0 | pv > 1, na.rm = TRUE)) {
      stop(sprintf("Column '%s' has values outside [0, 1]; if it holds -log10(P), pass it as `log10p =`.", col$p),
           call. = FALSE)
    }
    out$logp <- -log10(pv)
    zero <- which(pv == 0)
    if (length(zero) && (is.character(raw) || is.factor(raw))) {
      out$logp[zero] <- .logp_from_text(raw[zero])
      zero <- zero[is.na(out$logp[zero])]
    }
    if (length(zero)) {
      z <- out$beta[zero] / out$se[zero]
      ok <- is.finite(z)
      out$logp[zero[ok]] <- .logp_from_z(z[ok])
      if (!all(ok)) {
        out$logp[zero[!ok]] <- -log10(.Machine$double.xmin)
        warning(sum(!ok), " P-values of exactly 0 could not be recovered from ",
                "beta/se and were set to 1e-308.", call. = FALSE)
      }
    }
  }

  for (k in keep) out[[k]] <- data[[k]]

  ok <- out$chr %in% .chr_levels & is.finite(out$pos) & is.finite(out$logp)
  if (any(!ok)) {
    if (!quiet) message("Dropped ", sum(!ok), " variants with missing or non-standard chromosome, position or P-value.")
    out <- out[ok, , drop = FALSE]
  }
  rownames(out) <- NULL
  used <- unlist(col[!vapply(col, is.null, logical(1))])
  if (!quiet) message("Using columns: ", paste(sprintf("%s = %s", names(used), used), collapse = ", "))
  attr(out, "columns") <- used
  attr(out, "has_effect") <- any(is.finite(out$beta)) && any(is.finite(out$se))
  attr(out, "has_maf") <- any(is.finite(out$maf))
  attr(out, "effect_scale") <- if (is.null(col$beta) && !is.null(col$or)) "log_or" else "beta"
  attr(out, "effect_allele") <- if (effect_allele == "minor" && !is.null(col$eaf)) "minor" else "as_is"
  class(out) <- c("brownstone_sumstats", class(out))
  out
}
