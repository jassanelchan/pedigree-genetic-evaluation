# Project specification and execution

The JSON specification records the phenotype sources, model and software choice. Data paths may be absolute or relative to the specification's directory; executable paths should be absolute. Use forward slashes for Windows paths. The CLI runs in a fresh R session and does not require an existing R workspace.

## Top-level fields

| Field | Meaning |
| --- | --- |
| `schema_version` | Integer `1` |
| `engine` | Required explicit `blupf90`, `dmu`, or `compare` |
| `output_dir` | Output root; each new run gets a unique subdirectory |
| `software` | Overrides host profile: `renumf90`, `blupf90`, `dmu_dir` |
| `timeout` | Per-process maximum seconds; default 7200 |
| `threads` | Requested threads; default 2, bounded by host capacity |
| `pedigree` | Global pedigree source; may be overridden per analysis |
| `analyses` | Nonempty list of independently prepared models, with unique names |

`dmu_dir` points to the directory containing `dmu1.exe`, `dmuai.exe`, and `dmu4.exe`. R5.6 also requires its bundled shared runtime libraries. The skill does not redistribute these binaries. Read the backend guide before downloading from the official source. Package dependencies are `jsonlite`, `digest`, `processx`, `openxlsx`, `Matrix`, and `readxl` for Excel input.

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
| `model` | `animal`, `litter`, or `repeatability` |
| `traits` | Ordered list of one/two trait definitions |
| `fixed` | List of fixed class effects/covariates; an explicit mean is added |
| `shared` | For litter: `column` for verified litter ID, or `dam` plus `birth_date` |
| `event_date`, `birth_date` | Source column names/indices; optional unless used by effects or keys |
| `record_key` | Source columns uniquely identifying a record; required for repeated records |
| `qc` | `method:"sd"`, `sd:3`, optional `group_by` source columns; or explicit `method:"none"` |
| `inbreeding` | Boolean, default true |
| `predict_ids` | Optional original IDs without phenotypes that are represented in supplied pedigree |
| `allow_observed_founders` | Explicit authorization to add otherwise unrepresented observed/target IDs as unrelated founders; default false |
| `trend_basis` | `event_year` default, or explicit `birth_year` |
| `initial` | Optional named additive/shared/residual square nested arrays in trait order |
| `zero_covariances` | Explicit constraints `{component,trait1,trait2}`, indices one-based; component additive/shared/residual |

Traits need `name`, `column`, `unit`, optional `label` and `distribution`. Supported distributions are `gaussian` and `count_gaussian`; the latter validates nonnegative integers, then fits a Gaussian model. It is not a Poisson model. Unit/trait choice must come from user/data evidence; this pipeline does not manufacture weight corrections.

Fixed effects have `name`, `column`, `type` (`class` default or `cov`). Date-derived class effects instead use `from_date:"event"` or `"birth"` and `grouping:"year"`, `"month"`, or `"quarter"`. The same effects apply to both traits. Fixed-effect coverage is counted separately per trait, including levels with no observation for one trait. Rank checks use a reference-coded estimable basis; an empty cell alone is not a reason to exclude an entire year.

Use record keys including animal and event. Date fields used in keys are normalized before duplicate checks. Different date spellings cannot make one farrowing into two records. Repeated litter size should refer to the sow's phenotype, not the piglets' individual IDs. A permanent environment effect is different from the sow's own birth litter.

QC is a single pass on finite observations in the stated grouping. Boundaries are estimated before rejection, including sample SD with denominator n−1. Rejected trait values become missing; retain any record with at least one usable trait. Rows with no phenotype do not impose fixed-effect requirements, but supplied parents are preserved. Class and litter codes are rebuilt after filtering. Required fixed effects or litter IDs missing on usable records cause a targeted error; they are not silently imputed.

Default nonzero free covariance starts have correlation 0.1 and component fractions additive 0.3, shared 0.1 when present, residual 0.6/0.7. They are starting values, not final estimates. A zero covariance is a model constraint and must be declared separately. Version 1 supports off-diagonal zero constraints only; variance fixing is outside the generic interface.

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
