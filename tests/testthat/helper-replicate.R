# A second GWAS of the same variants (same IDs and alleles): every signal
# shrunk by `shrink`, with its own noise, so the null z stays N(0, 1).
replicate_gwas <- function(gw, shrink = 0.75, seed = 7) {
  withr_seed <- function(code) {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv())
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()))
    set.seed(seed)
    code
  }
  withr_seed({
    z <- shrink * gw$BETA / gw$SE + sqrt(1 - shrink^2) * stats::rnorm(nrow(gw))
    out <- gw
    out$BETA <- z * gw$SE
    out$P <- 2 * stats::pnorm(-abs(z))
    out
  })
}
