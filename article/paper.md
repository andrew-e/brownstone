---
title: 'maffect: A modern GWAS visualisation library for modern GWAS data'
tags:
  - R
  - genetics
  - genome-wide association study
  - GWAS
  - data visualization
  - Manhattan plot
authors:
  - name: Andrew Elmore
    orcid: 0000-0003-0292-447X
    affiliation: 1
affiliations:
  - name: University of Bristol, United Kingdom
    index: 1
    ror: 0524sp257
date: 28 October 2026
bibliography: paper.bib
---

# Summary

Genome-wide association studies (GWAS) test millions of genetic variants for association with a trait, and are usually summarised in a Manhattan plot of −log~10~ *P*-values by genomic position. That suits a few loci. Modern GWAS, biobank-scale and sequencing-based, include rare variants and find hundreds of significant loci, so *P* alone no longer says which matter. Readers need each signal's effect size, allele frequency and gene, and the *P*-value conflates the first two.

The R package maffect is a suite of plots for such data, built around a Manhattan plot that carries all three (\autoref{fig:sim}). Beneath the chromosome axis, an aligned track shows each signal's effect size with its 95% confidence interval. Significant variants are coloured by minor allele frequency (MAF), and loci are labelled with their nearest protein-coding gene, favouring the largest effects. Companion plots share the encoding (\autoref{fig:suite}): effect size against frequency, with the study's detection limit; a volcano plot; a regional plot with linkage disequilibrium (LD); Miami and forest plots for two GWAS; and a QQ plot by frequency class. The package reads common association tools' output, plots tens of millions of variants in seconds (\autoref{fig:finngen}), and flags lead variants whose effects disagree with their own case and control frequencies.

![A simulated GWAS of one million variants, from `manhattan(gwas, suggestive = 1e-5)`. Top: −log~10~ *P*, log-compressed above 20; solid points reach *P* < 5 × 10^−8^, hollow rings 10^−5^. Colour is the MAF class. Bottom: \|β\| of every hit, leads with 95% confidence intervals.\label{fig:sim}](figures/figure1_simulated.png){width="100%"}

# Statement of need

Effect size and allele frequency are central to interpreting GWAS results [@manolio2009missing]. Power depends on both, so rare variants reach significance only with large effects and common ones with small effects; a *P*-value alone conflates the two. Sequencing-based GWAS now include rare variants: FinnGen release 12 has 21.3 million per endpoint [@kurki2023finngen]. Readers need to see in one figure whether a signal is common or rare and how large; instead it is reported in tables or in separate effect-size-against-frequency ("trumpet") plots [@corte2023trumpet], which discard genomic position. maffect, for those who report or review GWAS results, keeps the familiar Manhattan layout and adds effect size, frequency and gene context on the same genome axis, from one call.

# State of the field

qqman [@turner2018qqman] made Manhattan plots routine in R. CMplot (in rMVP [@yin2021rmvp]), ggmanh [@lee2021ggmanh] and topr [@juliusdottir2023topr] add styling, axis rescaling, multi-trait views and gene annotation; topr also draws regional plots. LocusZoom [@pruim2010locuszoom] is the standard for single regions, and hudson [@lucas_hudson] draws mirrored ("Miami") plots. None places effect sizes or allele frequencies on the genome-wide figure. Manhattan++ [@grace2019manhattanpp] comes closest, with a transposed Manhattan heatmap beside a table of lead variants and their odds ratios and frequencies. Trumpet plots, in TrumpetPlots and GWASLab [@corte2023trumpet; @he2023gwaslab], show effect against frequency but not position.

We built a new package because an effect track sharing the genome axis needs an effect estimate, a standard error and a frequency for every point, not only a *P*-value, changing both layout and data model. It builds on ggplot2, patchwork and ggrepel [@wickham2016ggplot2; @pedersen_patchwork; @slowikowski_ggrepel], and returns editable plot objects.

# Software design

**Input.** `standardise_sumstats()` recognises the column names of PLINK, REGENIE [@mbatchou2021regenie], BOLT-LMM, SAIGE [@zhou2018saige], METAL and GWAS-SSF; any column can also be named. *P*-values below double precision are converted exactly on the log scale, and those that underflowed to zero are recovered from β/SE.

**Loci and labels.** Significant variants are clumped into loci by distance (±500 kb) and labelled with the nearest protein-coding gene (bundled Ensembl annotations [@dyer2025ensembl]) or the user's labels. When labels are limited, plots label alternately the most significant remaining locus and the one with the largest remaining effect, so elevated rare variants are named.

**Effect track.** The track shows \|β\|, the effect per copy of the trait-increasing allele (for odds ratios, of the risk allele), because the sign of one variant's effect only reflects which allele the GWAS reported; signed effects are optional. Leads are large points with confidence intervals. Allele frequencies bound LD: for minor allele frequencies *f*~1~ ≤ *f*~2~, *r* ≤ √(*f*~1~(1 − *f*~2~)/(*f*~2~(1 − *f*~1~))) [@vanliere2008r2]. If the lead is the only causal variant, a proxy's expected effect cannot exceed the lead's, so a point clearly above its lead suggests another signal. The same bound adds, without an LD panel, the strongest variant of each other frequency class in a locus when LD cannot explain it. Effects are not corrected for winner's curse [@kraft2008curses].

**Frequency, colour and axes.** Significant variants are coloured in four MAF classes (< 0.1%, 0.1–1%, 1–5%, ≥ 5%). Under simulated protanopia, deuteranopia and tritanopia [@machado2009cvd], every pair stays distinguishable (smallest ΔE = 9.2, OKLab × 100). The −log~10~ *P* axis is linear up to 2.5 times the significance line and log-compressed above it if the data go higher. The effect axis, in every plot that has one, is linear unless an estimate exceeds twice a break set at three times the median lead effect; it then breaks likewise, and the compressed zone always carries a tick.

**Effect against frequency.** `effect_vs_maf()` plots |β| against MAF on a log scale, with significance shown by style (solid, hollow or grey cloud). The curve is the smallest effect the study could detect at each frequency, with a chosen power (80% by default) at significance level α: (*z*~1−α/2~ + *z*~power~) × SE [@sham2014power], with typical standard error *c*/√(2*f*(1 − *f*)). The constant *c* absorbs the sample size, case–control balance and test calibration, so it is estimated as the median of SE × √(2*f*(1 − *f*)); no sample size is needed. An empty corner means "too small to see", not "not there". In FinnGen, *c* = 0.0053, and the smallest detectable odds ratio is 2.1 at MAF 0.1%, 1.27 at 1% and 1.05 at 50%.

**Artefact screen.** Some case–control summary statistics (FinnGen, SAIGE) report allele frequencies in cases and controls, which imply a crude odds ratio; covariate adjustment cannot move the reported estimate by an order of magnitude. `find_artifacts()` compares the two at each locus's lead variant and flags a disagreement of at least ten-fold or six standard errors; the plots ring flagged variants in black. It is a screen, not a test.

**Other plots.** `region_plot()` colours points by r² with the lead variant, from LDlink [@machiela2015ldlink] or the user's own table, and shapes them by frequency class (circle, square, triangle, star): a rare variant cannot be in strong LD with a common lead, so a pale triangle or star standing high is a signal of its own. The r² bins are one hue changing in lightness (smallest ΔE 13.3 under any simulated deficiency); the default bins of LocusZoom.js and locuszoomr fall to ΔE 5.3 and 4.0. `miami_plot()` draws two GWAS above and below the chromosome axis, and `forest_plot(compare = )` shows their signed effects side by side. `qq_plot()` draws one curve per MAF class and prints each genomic inflation factor, because rare-variant tests are often calibrated differently (in FinnGen, 1.25 for common variants and 1.05 below 0.1% MAF).

**Scale and quality.** Variants below the lowest threshold are thinned to one point per grid cell. On a 10-core Apple M4 laptop with 32 GB of memory, `manhattan()` built the plot of the 21.3 million FinnGen variants (\autoref{fig:finngen}) in 4.2 s, 1.8 s of which was column standardisation, and saving it as a PNG took 1.2 s. Reading the file took 18 s, and the run peaked at 8.9 GB, including the 3.15 GB input table. The package has 541 passing test expectations and passes `R CMD check` cleanly; a GitHub Actions workflow is set up to run it on Linux, macOS and Windows.

![FinnGen release 12, pure hypercholesterolaemia (E4_HYPERCHOL; 42,985 cases, 420,066 controls, 21.3 million variants), default settings. REGENIE log odds ratios are shown as \|β\|; the effect axis breaks at 0.5. Rings mark the two leads whose effects disagree with their own allele frequencies (near *PODN* and *MAN2B1*).\label{fig:finngen}](figures/figure2_finngen_hyperchol.png){width="100%"}

![The same data in `effect_vs_maf()`. The curve is the smallest effect detectable with 80% power at *P* < 5 × 10^−8^. Hits near it at the rare end passed largely because their estimates came out high.\label{fig:evm}](figures/figure3_finngen_effect_vs_maf.png){width="70%"}

![Other plots, on simulated data. A: `volcano_plot()`. B: `region_plot()` of *LDLR*, whose lead is a rare variant (triangle). C: `miami_plot()` of two separately simulated cohorts.\label{fig:suite}](figures/figure4_suite.png){width="80%"}

# Research impact statement

The package is new and not yet used in published work. Its significance rests on these tests, reproducible with the repository's scripts:

- **FinnGen release 12 pure hypercholesterolaemia** (\autoref{fig:finngen}): 72 loci at *P* < 5 × 10^−8^, with 8 lead variants at MAF < 0.1% and 6 at 0.1–1%.
- **GLGC 2013 LDL cholesterol** [@willer2013lipids]: 2.4 million variants, including *P* = 1.2 × 10^−652^ for rs7412.
- **Example outputs shipped with REGENIE and SAIGE** (binary traits, including SAIGE on chromosome X): the package's −log~10~ *P*-values matched the tools' exactly.

The plots also made two implausible estimates in the FinnGen data visible (\autoref{fig:evm}), and `find_artifacts()` flags exactly these two among the 72 leads. An ultra-rare variant near *MAN2B1* (MAF 0.009%, *P* = 5.8 × 10^−15^) has β = 12.6 (SE 1.61), an odds ratio of about 2.9 × 10^5^, while its case and control frequencies imply 1.75: a factor of 1.7 × 10^5^ (7.5 standard errors). An ultra-rare variant near *PODN* (MAF 0.04%, *P* = 7.3 × 10^−9^) has β = −4.76 (SE 0.82) against a crude log odds ratio of −0.30, a factor of 86 (5.4 standard errors). We did not examine why they are inflated; a conventional Manhattan plot shows them as moderate points. We expect it to be most useful for reporting sequencing-based GWAS and checking summary statistics.

# AI usage disclosure

The software, tests and documentation were produced with the assistance of a large language model (Anthropic Claude Opus 5.5, used through Claude Code), at the direction of the author. The author specified the goals, reviewed the figures at each design iteration, and reviewed and edited the code. 

Correctness was checked by unit tests and `R CMD check`; by comparing −log~10~ *P*-values with REGENIE and SAIGE outputs; by runs on public data (FinnGen, GLGC) with inspection of every figure; by checking the palettes against a published model of colour vision deficiency; and by verifying every reference against Crossref, DataCite, PubMed or the publisher's record.

The author reviewed and edited all generated code and text and takes responsibility for the content.

# Acknowledgements

We acknowledge the participants and investigators of the FinnGen study, and thank the Global Lipids Genetics Consortium for making their summary statistics public. [Funding information needed.]

# References
