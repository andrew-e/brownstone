gw <- simulate_gwas(n_variants = 5e4, seed = 1)

test_that("one family name is registered with every kind of device", {
  expect_equal(.family(), "maffect-sans")
  expect_true("maffect-sans" %in% systemfonts::registry_fonts()$family)
  expect_false(is.null(grDevices::pdfFonts("maffect-sans")[[1]]))
  expect_false(is.null(grDevices::postscriptFonts("maffect-sans")[[1]]))
  expect_null(.find_sans("No Such Font 123"))
  installed <- systemfonts::system_fonts()$family[1]
  expect_equal(.find_sans(c("No Such Font 123", installed)), installed)
})

test_that("every text the plots draw is in the package's family, labels included", {
  text_families <- function(plot) {
    unlist(lapply(plot$layers, function(l) {
      if (inherits(l$geom, c("GeomText", "GeomTextRepel"))) l$aes_params$family
    }))
  }
  b <- manhattan(gw, suggestive = 1e-5)
  expect_equal(b[[1]]$theme$text$family, "maffect-sans")
  expect_equal(unique(text_families(b[[1]])), "maffect-sans")
  expect_equal(unique(text_families(qq_plot(gw))), "maffect-sans")
  expect_equal(unique(text_families(volcano_plot(gw))), "maffect-sans")
  sim_ld <- attr(gw, "ld")
  causal <- attr(gw, "causal")
  h <- causal[causal$gene == "HMGCR", ]
  r <- suppressMessages(region_plot(gw, snp = paste0(h$chr, ":", h$pos), ld = sim_ld))
  expect_equal(unique(text_families(r[[3]])), "maffect-sans")
})

test_that("the family can be changed everywhere with one option", {
  op <- options(maffect.family = "")
  on.exit(options(op), add = TRUE)
  expect_equal(theme_maffect()$text$family, "")
  b <- manhattan(gw)
  expect_equal(b[[1]]$theme$text$family, "")
  repel <- Filter(function(l) inherits(l$geom, "GeomTextRepel"), b[[1]]$layers)[[1]]
  expect_equal(repel$aes_params$family, "")
})
