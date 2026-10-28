[English](provenance.md) | [简体中文](provenance.zh-CN.md)

# Synthetic data

All phenotype and pedigree files in this directory are simulated. They contain no farm data.

## Original shortcuts

[`project.json`](project.json) retains the bivariate growth and single-trait repeatability examples. Its pedigree includes founders, related descendants and a known inbreeding case. Phenotypes combine additive, shared and residual terms; total number born is additionally rounded to a nonnegative count.

## Mixed single and repeated records

[`multitrait-project.json`](multitrait-project.json) uses `mixed-records.csv` and `mixed-pedigree.csv` for a joint Gaussian model of age, backfat and total number born. The original files contain 1,776 event records and 780 pedigree rows: 720 animals each have one growth record, and 352 sows additionally have three farrowing records. Only sows have TNB observations; some growth records lack backfat. Growth values are never copied into farrowing events.

The model applies birth litter to the two growth traits, sow permanent environment to TNB, and different fixed effects to their respective traits. Genetic covariances span all three traits. Growth–TNB residual covariances are explicitly fixed to zero because those traits have no same-event observation. The specification uses `qc.method:"none"` to keep this synthetic example focused on model structure, not to recommend skipping QC on real data.

[`generate-general.R`](generate-general.R) reproduces the mixed files and specification with seed 20261028. Run it with `Rscript --vanilla` and optionally pass an output directory. The selected prediction pedigree is the target-and-ancestor closure, not necessarily all 780 original rows; unused founders are excluded.

Classification effects and event dates in both examples illustrate fixed effects and annual summaries. They do not represent a production population or recommended physiological ranges. The simulated component values are not reference estimates for a real breed.

Examples require a user-provided BLUPF90 or DMU installation. Executables are not distributed with the repository.
