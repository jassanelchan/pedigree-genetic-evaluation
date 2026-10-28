[English](project-spec.md) | [简体中文](project-spec.zh-CN.md)

# Project specification and execution

The JSON specification records the phenotype sources, model and software choice. Data paths may be absolute or relative to the specification's directory; executable paths should be absolute. Use forward slashes for Windows paths. The CLI runs in a fresh R session and does not require an existing R workspace.

## Top-level fields

| Field | Meaning |
| --- | --- |
| `schema_version` | Integer `1` or `2`; existing version-1 projects remain valid |
| `engine` | Required explicit `blupf90`, `dmu`, or `compare` |
| `output_dir` | Output root; each new run gets a unique subdirectory |
| `software` | Overrides host profile: `renumf90`, `blupf90`, `dmu_dir` |
| `timeout` | Maximum seconds per F90 process, or DMU1 plus its following module; default 7200 |
| `threads` | Requested threads; default 2 |
| `pedigree` | Global pedigree source; may be overridden per analysis |
| `analyses` | Nonempty list of independently prepared models, with unique names |

`dmu_dir` points to the directory containing `dmu1.exe`, `dmuai.exe`, and `dmu4.exe`. R5.6 also requires its bundled shared runtime libraries. The skill does not redistribute these binaries. Read the [backend guide](dmu.md) before downloading from the official source. Package dependencies are `jsonlite`, `digest`, `processx`, `openxlsx`, `Matrix`, and `readxl` for Excel input.

## Table source

```json
{
  "path": "phenotypes.xlsx",
  "sheet": "Sheet1",
  "columns": {"id": "个体号", "sire": "父亲", "dam": "母亲"},
  "missing_values": [-999],
  "date_origin": "1899-12-30"
}
```

- Only `columns.id` is required for phenotype sources. Optional sire/dam fields supplement real parent information for animals absent from the pedigree file. They are checked for contradictions.
- Pedigree sources require `columns.id`, `columns.sire`, and `columns.dam`.
- Optional `pedigree.supplements` is a list of additional table sources with verified ID/sire/dam mappings, such as real parent fields in a growth export. These sources supplement the pedigree for every analysis using that pedigree source; contradictory known parents fail. Declare shared-population sources explicitly rather than leaving an ancestor with known parents as a founder in another model.
- A column may be its exact name or a one-based integer position. Prefer exact names for Excel.
- Multi-sheet workbooks require a sheet choice. All fields are read as text to preserve letter IDs and leading zeros; previously lost Excel leading zeros cannot be recovered without another source.
- Text sources support `header` (default true), `sep` (CSV comma, TSV tab, otherwise whitespace), and `encoding` (default UTF-8-BOM). Use `header:false` for headerless native input files.
- Empty/NA/N/A/NULL are missing. Numeric trait values additionally recognize NAN and explicit `missing_values`. Parent IDs `0`, `-`, `.` are unknown.
- Numeric XLSX dates use verified workbook date-system metadata if `date_origin` is absent. Numeric dates in CSV/XLS need an explicit verified origin. Text dates accept YYYY-MM-DD or YYYY/MM/DD by default, or `date_format` using R date directives.

## Analysis fields

| Field | Meaning |
| --- | --- |
| `name`, `species` | Required analysis name and known species |
| `phenotypes` | Table source |
| `model` | `animal`, `litter`, `repeatability`, or `general` |
| `traits` | Ordered, nonempty list of trait definitions; no fixed two-trait limit |
| `fixed` | Fixed classes/covariates, optionally restricted to named traits; an explicit mean is added |
| `shared` | Litter shortcut only: `column` for verified litter ID, or `dam` plus `birth_date` |
| `random` | General model only: named independent `litter`, `pe` or `iid` terms |
| `event_date`, `birth_date` | Source column names/indices; optional unless used by effects or keys |
| `record_key` | Source columns uniquely identifying a record; required for repeated records |
| `qc` | `method:"sd"`, `sd:3`, optional `group_by` source columns; or explicit `method:"none"` |
| `inbreeding` | Boolean, default true |
| `predict_ids` | Optional original IDs without phenotypes that are represented in supplied pedigree |
| `allow_observed_founders` | Explicit authorization to add otherwise unrepresented observed/target IDs as unrelated founders; default false |
| `trend_basis` | `event_year` default, or explicit `birth_year` |
| `initial` | Optional full trait-by-trait nested arrays keyed by additive, random-term names and residual |
| `zero_covariances` | Explicit constraints `{component,trait1,trait2}`; component key and global one-based trait indices |

Traits need `name`, `column`, `unit`, optional `label`, `distribution`, `record_type` and `record_filter`. Supported distributions are `gaussian` and `count_gaussian`; the latter validates nonnegative integers, then fits a Gaussian model. It is not a Poisson model. Unit/trait choice must come from user/data evidence; this pipeline does not manufacture weight corrections.

Fixed effects have `name`, `column`, `type` (`class` default or `cov`). Date-derived class effects instead use `from_date:"event"` or `"birth"` and `grouping:"year"`, `"month"`, or `"quarter"`. Optional `traits:["age120","bf120"]` restricts an effect to those trait names; omission means all traits. For example, sex can apply to growth traits while parity applies only to litter size. Fixed-effect coverage and rank are checked separately for each applicable trait. A fixed value is required only on rows with an observed applicable trait; inactive rows use neutral placeholders internally, not imputed observations. Rank checks use a reference-coded estimable basis; an empty cell alone is not a reason to exclude an entire year.

## Model shortcuts and named random terms

All models have one pedigree-based additive genetic term on every trait and one residual covariance matrix. The existing shortcuts accept one or more traits:

- `animal`: no additional random term.
- `litter`: one identity-related common-litter term named `shared`, defined through `shared`.
- `repeatability`: one animal permanent-environment term named `shared`; traits default to repeated records.

Use `general` to combine independent terms on different trait subsets. Do not include `shared` in this model. An example random-effects section is:

```json
"random": [
  {"name": "birth_litter", "role": "litter", "column": "litter", "traits": ["age120", "bf120"]},
  {"name": "sow_pe", "role": "pe", "traits": ["TNB"]}
]
```

The additive term is implicit; do not add it to `random`. Each additional term needs a unique name matching `[A-Za-z][A-Za-z0-9_]*`; `additive`, `residual` and `shared` are reserved. `traits` is optional and defaults to all traits. Roles determine grouping and interpretation:

| Role | Grouping | Interpretation |
| --- | --- | --- |
| `litter` | Verified grouping `column`, or `dam` with analysis `birth_date` | Common birth-litter effect |
| `pe` | Animal ID; no alternative `column` or `dam` allowed | Permanent environment across repeated observations |
| `iid` | Required grouping `column` | Other identity-related grouped random effect |

Each term can have covariance among its active traits; different named terms are independent blocks. Common litter and generic `iid` variance are not included in the repeatability numerator merely because they are environmental. Each identity-related term needs at least two observed groups and at least one replicated group for every active trait. Identical grouping partitions for two terms on the same observed trait are rejected; additional terms still require a scientifically justified, identifiable model.

Current scope is Gaussian pedigree models. Genomic relationships, maternal genetic effects, unknown-parent genetic groups, heterogeneous residuals, threshold models and random regression are not implemented. More traits and covariance parameters also require more data and computation; a larger model is not automatically better.

## Single and repeated records in one model

`record_type` is `single` by default, or `repeated`; the repeatability shortcut defaults to `repeated`. A single-record trait may have at most one observed value per animal. Any repeated trait requires a `record_key` containing the animal ID and real event identifiers or dates. Never copy a growth phenotype into each farrowing row to create repeated records.

For a general model mixing single and repeated traits, every trait must declare its event membership. For example:

```json
{"name": "age120", "column": "age120", "unit": "days", "record_type": "single", "record_filter": {"column": "event_type", "values": ["growth"]}}
```

Use the corresponding `farrowing` filter for a repeated TNB trait. A filter defines which rows belong to a trait, including rows where that trait was not measured. An observed value outside its declared filter is an error, not a value silently discarded. Keep traits that do not belong to a row missing, and leave inapplicable fixed/random-effect fields empty where appropriate.

The [multitrait example](example/multitrait-project.json) combines two single-record growth traits and repeated TNB. Growth and farrowing are separate real event rows; sex applies to growth, parity to TNB, birth litter to growth, and sow PE to TNB. Genetic covariances connect the traits through the animals and pedigree. No row observes growth and TNB together, so their residual covariances are explicitly fixed to zero in `zero_covariances`.

Use record keys including animal and event. Date fields used in keys are normalized before duplicate checks. Different date spellings cannot make one farrowing into two records. Repeated litter size should refer to the sow's phenotype, not the piglets' individual IDs. A permanent environment effect is different from the sow's own birth litter.

QC is a single pass on finite observations in the stated grouping. Boundaries are estimated before rejection, including sample SD with denominator n−1. Rejected trait values become missing; retain any record with at least one usable trait. Rows with no phenotype do not impose fixed-effect requirements, but supplied parents are preserved. Class and random-group codes are rebuilt after filtering. Required fixed effects or random groups missing on applicable observed rows cause a targeted error; they are not silently imputed.

## Covariance matrices and constraints

Every component matrix in `initial` uses the complete ordered trait list, even when a term applies to a subset. Keys are `additive`, `residual` and the declared random-term names; shortcuts retain the name `shared`. With traits `[age120,bf120,TNB]`, a litter term active on the first two traits has a 3×3 matrix whose last row and column are zero. The active principal submatrix must be positive definite. Inactive entries are structural zeros and are excluded from estimation.

For user constraints, `trait1` and `trait2` must be integers referring to the global one-based trait order, not a backend's local numbering. For example:

```json
"zero_covariances": [
  {"component": "residual", "trait1": 1, "trait2": 3},
  {"component": "residual", "trait1": 2, "trait2": 3}
]
```

When no retained real record observes a trait pair together, their residual covariance is not estimable and must have an explicit, justified zero constraint. An identity-related term likewise requires a zero covariance for active trait pairs that have no observed group in common, even across separate records. The workflow stops rather than choosing these assumptions automatically. Genetic covariance may still be estimated through connected animals and pedigrees, subject to data support. Declaring a residual zero does not fix the genetic covariance to zero.

Default starts allocate 0.3 of each trait's sample variance to additive genetics. Where additional terms apply, they divide a total fraction of 0.1, leaving 0.6 for residuals; without them the residual fraction is 0.7. A block with k active traits starts with free correlations `0.1/max(1,k−1)`. These are numerical starting values, not final estimates. Free covariance starts must be nonzero; fixed-zero covariances must be declared separately. Only off-diagonal zero constraints are supported; fixing an active variance is outside the generic interface.

## Commands

```powershell
$env:LC_ALL = 'English_United States.utf8'
$env:LANG = 'English_United States.utf8'
& '<Rscript-path>' --vanilla '<skill-path>\scripts\run_eval.R' --spec '<project-path>\project.json'
```

Use the discovered Rscript path. `--engine` and `--output` override their specification fields. `--prepare-only` prepares without running either solver. `--resume <prior-run-directory>` reads verified caches only when fingerprints match. Scripts and the current specification are copied into a new invocation directory on every call. Data/model fingerprints deliberately exclude trend-only year grouping so changing graphs does not trigger REML.

## Output and failure state

`delivery/evaluation_summary.xlsx` combines models, descriptive statistics/CV, QC, fixed-level coverage, input audit, variance components/parameters, annual means and slopes. Per-engine folders contain original-ID EBV CSV, detailed workbook and figures. Comparison folders include component, EBV/PEV and fitted-value comparisons where native output permits.

`analysis_<n>/` preserves prepared context, ID map, rank checks, exact model definition and source hashes. Backend folders preserve all generated cards/data, REML attempts, logs, process status and verified fit/EBV caches. `code_records/` holds readable complete source snapshots. `run_manifest.json` records software/R/package provenance and states `prepared`, `complete`, or `partial`. Exit status 2 means inspect errors/partial output; do not announce a successful complete evaluation.

`input_sources/` preserves SHA256-verified copies of every raw input and parent supplement; unchanged inputs are deduplicated by hash. `input_manifest.csv` records original and snapshot paths. `reproduce_project.json` uses relative snapshot paths, so a transferred run does not depend on the original Desktop folder. Use a saved `code_records/invocation_.../run_eval.R` with that specification, and configure valid software paths/packages on the destination host. The pipeline never overwrites source files. Raw snapshots contain private herd data and belong to the project reproduction directory, never the shareable skill package.
