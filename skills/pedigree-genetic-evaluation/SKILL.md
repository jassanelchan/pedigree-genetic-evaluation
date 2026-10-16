---
name: pedigree-genetic-evaluation
description: Run pedigree-based livestock genetic evaluation with BLUPF90, DMU, or both. Supports one/two-trait Gaussian animal and common-litter models, single-trait repeatability, QC, variance components, EBV and annual trends. Use when evaluation is requested; method questions alone do not start computation. Excludes genomic, threshold and random-regression models.
---

# Pedigree genetic evaluation

**English** | [简体中文](SKILL.zh-CN.md)

Evaluate supplied phenotype and pedigree files, keeping reproducible inputs, code, parameter files and logs. Follow the software, trait definitions, model and QC rules agreed for the current analysis.

## Establish the analysis

For a question about methods or feasibility, explain without running data. For a planning request, provide the plan. Start the evaluation when analysis is requested.

Identify the phenotype and pedigree sources, species, traits, units and engine (`blupf90`, `dmu` or `compare`). Inspect sheets, column names and a few records together. Resolve file formats and column mappings from the files and existing requirements where possible.

Clarify ambiguities that change the analysis: uncertain traits or units, mixed populations, unclear date sources, conflicting parentage, duplicate events, unreliable litter definitions, or a proposed change to QC or model. Carry forward decisions already confirmed.

Use `pedigree.supplements` to combine genuine parent information from additional files. Preserve source records and checksums. Do not discard known relationships merely because one export lacks the parent rows.

## Configure the model

Read the [project configuration guide](references/project-spec.md) and write a JSON specification. The [synthetic example](references/example/project.json) provides a working structure.

- `animal`: one or two traits, with fixed effects, additive genetic effects and residuals.
- `litter`: one or two traits, adding common litter. Use a verified litter ID, or dam ID together with offspring birth date; dam ID alone does not distinguish litters.
- `repeatability`: one trait with repeated records, adding permanent environment linked by animal ID. Require a record key composed of animal and event.

Both traits use the same fixed-effect structure. Categorical effects and continuous covariates are supported; inbreeding is included by default. Genomic relationships, three-trait models, unknown-parent genetic groups, maternal genetic effects, heterogeneous residuals, threshold models and random regression are outside the supported scope. Explain a scope mismatch rather than silently dropping a requested component.

Interpret random effects by their role in the model. A software heading such as “Genetic variance for effect 4” does not establish that the effect is additive genetic. Do not hard-code the additive effect number when parsing results.

## Run and verify

Read the [BLUPF90 guide](references/blupf90.md) or [DMU guide](references/dmu.md) for the chosen engine; read both for a comparison.

1. Check R dependencies, executable paths and versions. Read paths from the project or local `host_profile.json`; keep personal paths out of public files. Do not substitute another engine when the selected program is unavailable.
2. Launch `scripts/run_eval.R` in a fresh process with `Rscript --vanilla`. Sources are UTF-8. If Windows requires an explicit locale, use `English_United States.utf8`. Do not load an old `.RData` or depend on scripts containing manually entered estimates.
3. Check IDs, parentage, dates, record keys and the fixed-effect design. Default QC is a single mean ±3 SD pass per trait, retaining records with at least one valid trait. Count missing and rejected values separately. Do not add a weight range or correction formula without an analysis requirement.
4. Run REML, read the final components, then use those components in BLUP to obtain EBV and SEP. Estimates must come from the actual program output.
5. Verify convergence, covariance matrices, parameter order, complete animal-by-trait coverage and original-ID mapping. After an initial numerical failure, make at most two further attempts with the same data and model. Reconfirm changes to QC, effect removal or model splitting.
6. Inspect `run_manifest.json`, the summary workbook, plots and logs before delivery. Keep an incomplete run marked partial when standard errors or required outputs are missing. Never replace the full sampling covariance matrix with a diagonal matrix made from squared component standard errors.

`--prepare-only` prepares and checks inputs without running solvers. `--resume` reuses verified fits, predictions and comparisons when fingerprints match. A change limited to plots or annual summaries does not require another variance-component fit.

## Summarise and deliver

Read the [result definitions](references/results.md). Keep full numerical precision internally; display phenotype statistics and CV to two decimals, counts as integers and genetic parameters to three decimals. Calculate parameter point estimates from the final components. Use the full sampling covariance matrix and the Delta method for their standard errors, distinguishing these from native random-sampling estimates.

In `compare` mode, use the same data and model. Check component estimates first, then compare EBV, PEV and fitted values with components held fixed at the same values. Account for recorded numerical and printing precision. A high correlation does not replace a difference check; raw likelihood constants from different programs are not directly comparable.

Default annual summaries to the observation event: measurement or end-test year for growth, farrowing year for reproduction. Count each animal's EBV once per year. Use birth cohorts only when requested and state the grouping basis. Raw phenotype trends also reflect environment and population composition; they do not measure genetic improvement alone.

Deliver the result workbook, original-ID EBV tables, QC and trend plots, and the complete reproduction directory. Preserve checked source copies, a relative-path reproduction specification, code, parameter files, ID maps, logs and versions for every run. Do not include real datasets or third-party executables in the skill package. Create a presentation or PDF only when requested.

Explain estimates using the actual model and statistical definitions. Heritability is a population variance ratio, not the fraction of an individual animal's deviation attributable to genetics. Agreement with a published range alone does not validate the calculation.
