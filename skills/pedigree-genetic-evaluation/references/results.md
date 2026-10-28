[English](results.md) | [简体中文](results.zh-CN.md)

# Results, uncertainty and figures

The result module operates on validated post-QC records and final model components.
It does not select animals, change the model or repair failed estimation. Each
random term has a unique component name and a semantic role: `additive`, `litter`,
`pe` or `iid`; the residual is separate. For example, `sow_pe` is a component name
whose declared role is `pe`. Names and roles are not inferred from the software's
effect numbers or printed headings. The legacy shortcuts retain their original
result labels while recording the source component `shared`.

## Genetic parameters

For trait t, the reported model-scale phenotypic variance is

`Vp[t] = additive[t,t] + sum(random_term[t,t]) + residual[t,t]`,

where the sum includes every declared additional independent random term
applicable to t. Inactive rows/columns are structural zero. This is variance after
the specified fixed effects, not the raw sample variance. The usual population parameter uses
additive variance rather than multiplying every animal by its own `1 + F`.

- Narrow-sense heritability: `h2[t] = additive[t,t] / Vp[t]`.
- Repeatability, only for traits with permanent environment:
  `(additive[t,t] + sum(applicable PE variances)) / Vp[t]`.
  Only terms declared with role `pe` enter this numerator; litter and `iid`
  terms remain in the denominator.
- Additive genetic correlation for every distinct pair s,t:
  `additive[s,t] / sqrt(additive[s,s] * additive[t,t])`.
- All component variances and cross-trait covariances are included in the same table,
  with units and standard errors. Ratios and correlations are dimensionless.

For the [three-trait example](example/multitrait-project.json), the two growth
traits have Vp=VA+Vbirth_litter+VE, while TNB has Vp=VA+Vsow_pe+VE.
Repeatability is reported for TNB only. The full additive matrix links all three
traits; the litter matrix has no TNB variance, and the PE matrix has no growth
variance. The parameter table records both `component` and semantic `role` so
several named environmental terms remain distinguishable.

A zero total variance makes heritability and repeatability undefined. A zero
additive diagonal makes the corresponding genetic correlation undefined. Small or
boundary estimates are retained and labeled; they are not replaced with familiar
literature values. Non-finite, asymmetric or materially non-positive-semidefinite
final matrices fail validation rather than yielding a normal-looking report.

## Standard errors

Standard errors of functions use the analytic Delta method:

`SE(f) = sqrt(gradient(f)' %*% param_cov %*% gradient(f))`.

`param_meta` must describe every estimated covariance parameter in the exact row
and column order of the full sampling covariance matrix. Its `role` field contains
the unique component key, with semantic roles supplied by the model registry;
trait indices follow the global trait order. Native factor-local indices are
mapped before use, including nonconsecutive active subsets. Effect order is not
assumed. Estimated parameters may be rearranged if metadata and covariance are
rearranged together. Constrained and inactive structural zeros are absent from
the sampling parameter set.

For heritability A/P, derivatives are `(P-A)/P^2` for additive variance and
`-A/P^2` for each other diagonal variance. For repeatability N/P, the derivative
is `(P-N)/P^2` for additive and every applicable permanent-environment variance, and `-N/P^2`
for the other diagonals. For genetic correlation r, derivatives are
`1/sqrt(A_s*A_t)` for the covariance, `-r/(2*A_s)` and `-r/(2*A_t)` for the
two additive variances.

The full covariance, including off-diagonal terms, is necessary. Do not replace a
missing matrix with independent component SEs. Missing, incomplete, asymmetric or
indefinite sampling covariance leaves derived SE unavailable with a reason. Native
component SE matrices may still be retained when verified, but cannot supply the
missing covariance terms for heritability or correlation. Boundary Delta SEs are
asymptotic approximations and are labeled accordingly. These SEs may differ slightly
from software Monte Carlo SEs because the methods differ.

## Reliability

Reliability is optional and is calculated only from verified SEP, additive variance
and the diagonal of the animal relationship matrix represented by verified F:

`REL = 1 - SEP^2 / ((1 + F) * additive[t,t])`.

The SEP must be the prediction-error standard deviation, not the sampling SE of a
variance component. F must be verified against the fitted relationship convention.
Missing F is not silently set to zero. Negative SEP, invalid F, zero additive
variance, or PEV greater than the corresponding model variance leave REL missing
with a status. Invalid values are not clipped to the interval [0,1]. This conditional
reliability does not include uncertainty from estimating the variance components.

## Annual summaries

The default `event_year` groups by the observation/production event. Growth uses
the measurement year, and litter size uses the farrowing year. Optional `birth_year`
summaries answer a different question: change across birth cohorts. The two bases
are identified in the table and filename; they are never silently mixed.

For mixed single/repeated designs, `record_filter` defines each trait's event
membership. Only retained rows belonging to that trait enter its yearly groups,
including applicable rows where that trait is missing but another is observed.
A sow's farrowing dates therefore do not create new growth-year EBV entries.
Membership is distinct from whether an observation passed QC.

For each year and trait:

- Phenotypic means and SD use valid observations for that trait. Repeated records
  contribute as records, so litter size has one observation per retained litter.
- Mean EBV counts each animal at most once per year within its applicable event
  rows, including animals with another trait observed but their own value missing.
  The same sow may enter several farrowing years for TNB, once in each year.
- `n_records`, `n_animals`, `n_phenotyped` and `n_ebv` are separate. An EBV yearly
  group is not weighted by each sow's number of litters.
- Invalid or unavailable years are omitted from that annual summary and reported.
  Conflicting birth years for the same animal fail validation.

At least three effective annual means are required for a linear trend. The fit is
ordinary unweighted regression of annual means on year. Output includes slope,
slope SE, R-squared, effective year count and trait units per year. Two-year groups
retain their means without claiming a fitted trend. No fit changes any EBV or
re-estimates the model. Changing the annual basis only requires result processing.

## Deliverables and presentation

`ge_export_results` creates a per-analysis `results.xlsx`, `EBV.csv`, QC histograms
and phenotype/genetic trend PNGs. It returns paths and the parameter/annual/trend
tables for the combined project summary.

- Numerical cells retain full precision. Excel display formats use integers for
  counts, two decimals for descriptive statistics and three for components,
  heritability and other genetic parameters.
- `CV_percent` is `100 * SD / mean` for a positive mean. CV is undefined when the
  mean is zero or negative; do not report it as an ordinary positive percentage.
- CSV has a UTF-8 BOM. Original IDs are strings; the workbook preserves leading
  zeros. Import CSV ID columns as text if Excel automatically converts them.
- QC before/after histograms use the same bins and axes. Provide optional
  `ctx$raw_data` or `ctx$qc_before` with numeric `y1 ... yN` values before QC.
  If these values are unavailable, export retained distributions only and report
  that before-QC figures could not be produced.
- Trend plot titles are only `表型进展` and `遗传进展`. Axes identify the annual
  basis, trait and unit. Lines show means; dashed lines show the annual fit.
  No explanatory footer is added. More than four traits are split across PNGs
  in blocks of at most four, with the trait range included in each filename.
- No valid values for a trait produce an explicitly empty panel, not a fictitious
  distribution. Report labels follow supplied trait names and units.
- Optional user descriptive tables should expose counts as integer columns,
  and CV as a column clearly indicating percent; the exporter does not reinterpret
  a user's supplied fractional CV without its definition.

Large EBV outputs exceeding the worksheet row limit remain complete in CSV rather
than being truncated to fit Excel. Source input files and model outputs are never
overwritten by this module.
