---
name: pedigree-genetic-evaluation
description: Run pedigree-based livestock genetic evaluation with BLUPF90, DMU, or both. Supports multitrait Gaussian animal, litter, repeatability and general models with trait-specific effects and mixed single/repeated records, QC, genetic parameters, EBV and annual trends. Use when evaluation is requested; method questions alone do not start computation. Excludes genomic, threshold and random-regression models.
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

Read the [project configuration guide](references/project-spec.md) and write a JSON specification. Schema versions `1` and `2` are accepted. The [original synthetic example](references/example/project.json) retains separate growth and reproduction analyses; the [multitrait example](references/example/multitrait-project.json) combines single growth measurements and repeated litter-size records.

- `animal`: one or more traits, with fixed effects, additive genetic effects and residuals.
- `litter`: one or more traits, adding common litter. Use a verified litter ID, or dam ID together with offspring birth date; dam ID alone does not distinguish litters.
- `repeatability`: one or more repeated traits, adding permanent environment linked by animal ID.
- `general`: one or more traits, with trait-specific fixed effects and multiple named random terms. Declare each term's role (`litter`, `pe` or `iid`), grouping and applicable traits. Additive genetic effects are implicit and apply to all traits. Named random terms are independent of each other; their within-term covariances apply only to their declared traits.

The shortcuts preserve existing project files. Use `fixed[].traits` to limit a class effect or covariate to selected trait names; omission applies it to all traits. Use `general` when different traits require different random terms, or single and repeated records are combined. PE must be grouped by animal, with actual repeated observations for the traits to which it applies. Do not rename a litter or other grouping term as PE merely to obtain repeatability.

Declare `record_type` and, for every trait in a mixed single/repeated design, an explicit `record_filter` identifying its real event rows. Repeated designs need a record key containing animal ID and event identifiers. Never copy a single growth observation into multiple farrowing rows. A trait may be applicable but unmeasured on a row; that differs from a row outside its event filter. If two traits have no paired observations, their residual covariance requires an explicit, scientifically justified zero constraint; do not assume it automatically.

All covariance arrays use the full trait order. Inactive rows and columns of a named term are structural zeros, not failed estimates. Free covariance starts must be nonzero, while user-specified zero constraints must be explicit. Validate the active covariance blocks, fixed-effect design and grouping before fitting; increasing trait count alone does not establish identifiability.

Inbreeding is included by default. Genomic relationships, unknown-parent genetic groups, maternal genetic effects, covariance between different named random terms, heterogeneous residuals, threshold models and random regression are outside the supported scope. Explain a scope mismatch rather than silently dropping a requested component. Verify support against the chosen adapter and program version, not the software family's capabilities in general.

Interpret random effects by their role in the model. A software heading such as “Genetic variance for effect 4” does not establish that the effect is additive genetic. Do not hard-code the additive effect number when parsing results.

## Run and verify

Read the [BLUPF90 guide](references/blupf90.md) or [DMU guide](references/dmu.md) for the chosen engine; read both for a comparison.

1. Check R dependencies, executable paths and versions. Read paths from the project or local `host_profile.json`; keep personal paths out of public files. Do not substitute another engine when the selected program is unavailable.
2. Launch `scripts/run_eval.R` in a fresh process with `Rscript --vanilla`. Sources are UTF-8. If Windows requires an explicit locale, use `English_United States.utf8`. Do not load an old `.RData` or depend on scripts containing manually entered estimates.
3. Check IDs, parentage, dates, record keys and the fixed-effect design. Default QC is a single mean ±3 SD pass per trait, retaining records with at least one valid trait. Count missing and rejected values separately. Do not add a weight range or correction formula without an analysis requirement.
4. Run REML, read the final components, then use those components in BLUP to obtain EBV and SEP. Estimates must come from the actual program output.
5. Verify convergence, active covariance blocks, named-component and global-trait parameter mapping, complete animal-by-trait coverage and original-ID mapping. Never infer a term's role from its name or a native local trait number. After an initial numerical failure, make at most two further attempts with the same data and model. Reconfirm changes to QC, effect removal or model splitting.
6. Inspect `run_manifest.json`, the summary workbook, plots and logs before delivery. Keep an incomplete run marked partial when standard errors or required outputs are missing. Never replace the full sampling covariance matrix with a diagonal matrix made from squared component standard errors.

`--prepare-only` prepares and checks inputs without running solvers. `--resume` reuses verified fits, predictions and comparisons when fingerprints match. A change limited to plots or annual summaries does not require another variance-component fit.

## Summarise and deliver

Read the [result definitions](references/results.md). Keep full numerical precision internally; display phenotype statistics and CV to two decimals, counts as integers and genetic parameters to three decimals. Calculate parameter point estimates from the final components. Use the full sampling covariance matrix and the Delta method for their standard errors, distinguishing these from native random-sampling estimates.

Report each named random component separately. For each trait, phenotypic variance sums all applicable diagonal components and residual variance; narrow-sense heritability uses additive variance in the numerator. Report repeatability only for traits with PE, using additive plus applicable PE variances in its numerator. Do not include litter or other `iid` terms there. Report additive genetic correlations for trait pairs using the global trait order.

In `compare` mode, use the same data and model. Check component estimates first, then compare EBV, PEV and fitted values with components held fixed at the same values. Account for recorded numerical and printing precision. A high correlation does not replace a difference check; raw likelihood constants from different programs are not directly comparable.

Default annual summaries to the observation event: measurement or end-test year for growth, farrowing year for reproduction. Respect each trait's event membership and count each animal's EBV once per trait and year, including applicable rows with a missing phenotype. Do not use a sow's later farrowing years as years for its growth-trait trend. Use birth cohorts only when requested and state the grouping basis. Raw phenotype trends also reflect environment and population composition; they do not measure genetic improvement alone.

Deliver the result workbook, original-ID EBV tables, QC and trend plots, and the complete reproduction directory. Preserve checked source copies, a relative-path reproduction specification, code, parameter files, ID maps, logs and versions for every run. Do not include real datasets or third-party executables in the skill package. Create a presentation or PDF only when requested.

Explain estimates using the actual model and statistical definitions. Heritability is a population variance ratio, not the fraction of an individual animal's deviation attributable to genetics. Agreement with a published range alone does not validate the calculation.
