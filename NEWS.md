# maffect 0.1.0

First release.

* `manhattan()` draws a Manhattan plot with an aligned effect-size track
  (with 95% confidence intervals), variants coloured by minor allele frequency,
  nearest-gene labels (GRCh37 and GRCh38), an optional suggestive tier, and
  y axes that break only when a locus or an effect far exceeds the rest.
* Companion plots with the same encoding: `effect_vs_maf()` (with the study's
  detection limit), `volcano_plot()`, `region_plot()` (linkage disequilibrium
  from LDlink or a table of your own), `miami_plot()`, `forest_plot()` and
  `qq_plot()` (one curve per frequency class).
* `find_artifacts()` flags lead variants whose reported effect disagrees with
  the crude effect implied by the case and control allele frequencies.
* `standardise_sumstats()` recognises PLINK, REGENIE, BOLT-LMM, SAIGE, METAL,
  FinnGen and GWAS-SSF columns; `find_loci()`, `nearest_gene()` and
  `ld_proxy()` are exported for use on their own.
* `simulate_gwas()` simulates summary statistics for trying the plots out.
