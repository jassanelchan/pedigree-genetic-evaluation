# Pedigree genetic evaluation

**English** | [简体中文](README.zh-CN.md)

A livestock genetic evaluation workflow built around R, BLUPF90 and DMU. Start with phenotype and pedigree files; estimate variance components, predict breeding values and summarise annual trends.

Choose BLUPF90, DMU, or run both against the same data and model. Each run keeps the input files, analysis code, parameter files, ID maps and solver logs.

## Get started

Use the workflow as a skill, from the command line, or from R / RStudio. All three use the same R scripts and need the [software setup](#software-setup) below for a full evaluation.

### 1. Import the skill

With Node.js and npm installed, copy this command into a terminal:

```sh
npx skills add jassanelchan/pedigree-genetic-evaluation --skill pedigree-genetic-evaluation
```

Choose your application and installation location when prompted. The installer adds the skill instructions and supporting files; R and the evaluation programs must be installed separately.

Then provide the phenotype and pedigree files with an analysis request, for example:

> Use $pedigree-genetic-evaluation with BLUPF90. Fit a bivariate model for age at 120 kg and backfat, with sex and measurement year-month as fixed effects and additive genetic and common litter effects as random effects. Use one-pass 3 SD QC and retain the code, parameter files and logs.

See the [skill instructions](skills/pedigree-genetic-evaluation/SKILL.md) for supported models and checks. The import command uses the [skills installer](https://github.com/vercel-labs/skills).

### 2. Run from the command line

Clone the repository and enter its directory:

```sh
git clone https://github.com/jassanelchan/pedigree-genetic-evaluation.git
cd pedigree-genetic-evaluation
```

Install R, make `Rscript` available on `PATH`, and install the dependencies once:

```sh
Rscript --vanilla -e "install.packages(c('jsonlite', 'digest', 'processx', 'openxlsx', 'readxl', 'Matrix'), repos='https://cloud.r-project.org')"
```

Check the bundled synthetic example first. This prepares the inputs without starting either solver:

```sh
Rscript --vanilla skills/pedigree-genetic-evaluation/scripts/run_eval.R --spec skills/pedigree-genetic-evaluation/references/example/project.json --engine blupf90 --prepare-only
```

After configuring BLUPF90, run the evaluation:

```sh
Rscript --vanilla skills/pedigree-genetic-evaluation/scripts/run_eval.R --spec skills/pedigree-genetic-evaluation/references/example/project.json --engine blupf90
```

Replace `blupf90` with `dmu` or `compare` to use DMU or both programs. For your own data, copy the example JSON, edit its column mappings and model, and pass its path to `--spec`. See the [project configuration guide](skills/pedigree-genetic-evaluation/references/project-spec.md).

### 3. Run from R / RStudio

Set the working directory to the cloned repository. If needed, install the dependencies in the R console:

```r
install.packages(
  c("jsonlite", "digest", "processx", "openxlsx", "readxl", "Matrix"),
  repos = "https://cloud.r-project.org"
)
```

Then run the example:

```r
script <- normalizePath(
  "skills/pedigree-genetic-evaluation/scripts/run_eval.R"
)

project <- normalizePath(
  "skills/pedigree-genetic-evaluation/references/example/project.json"
)

rscript <- file.path(
  R.home("bin"),
  if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
)

prepare_only <- TRUE

status <- system2(
  rscript,
  c(
    "--vanilla", shQuote(script),
    "--spec", shQuote(project),
    "--engine", "blupf90",
    if (prepare_only) "--prepare-only"
  )
)

stopifnot(status == 0L)
```

After configuring the solver, set `prepare_only <- FALSE` to run the full evaluation. This launches a separate R process; do not use `source()` on the command-line entry script.

## Software setup

The verified environment is **Windows with R 4.5.3**. Other operating systems and program versions have not been validated.

| Program | Verified version | Required for |
| --- | --- | --- |
| RENUMF90 | 1.169 | BLUPF90 and comparison runs |
| BLUPF90+ | 2.76 | BLUPF90 and comparison runs |
| DMU | Version 6, release 5.6 | DMU and comparison runs |

Install the required programs separately: [BLUPF90](https://nce.ads.uga.edu/wiki/) and [DMU](https://dmu.ghpc.au.dk/dmu/). Executables are not included. The DMU adapter uses the release with separate `dmu1`, `dmuai` and `dmu4` programs.

In the cloned repository, copy the local configuration template:

```powershell
Copy-Item skills/pedigree-genetic-evaluation/host_profile.example.json skills/pedigree-genetic-evaluation/host_profile.json
```

Open `host_profile.json` in a text editor and replace the example paths. Keep the entries needed by your chosen engine:

```json
{
  "software": {
    "renumf90": "C:/software/BLUPF90/renumf90.exe",
    "blupf90": "C:/software/BLUPF90/blupf90+.exe",
    "dmu_dir": "C:/software/DMU/bin"
  }
}
```

For an imported skill, make this copy inside the installed skill directory instead. `host_profile.json` is excluded from Git. A project's `software` settings override the local configuration. RENUMF90 and BLUPF90+ can also be found through `PATH`.

If a Windows session misreads Chinese paths or columns, set these variables in PowerShell before starting R:

```powershell
$env:LC_ALL = 'English_United States.utf8'
$env:LANG = 'English_United States.utf8'
```

## Models

| Model | Traits | Random effects |
| --- | --- | --- |
| Animal (`animal`) | 1 or more | Additive genetic, residual |
| Common litter (`litter`) | 1 or more | Additive genetic, common litter, residual |
| Repeatability (`repeatability`) | 1 or more | Additive genetic, permanent environment, residual |
| General (`general`) | 1 or more | Additive genetic, multiple named litter, permanent-environment or independent grouping effects, residual |

The first three are shortcuts for common designs; existing project files remain valid. The general model also supports single and repeated traits in the same analysis. Schema versions `1` and `2` are accepted. The number of traits is subject to data, identifiability and backend checks; three- and four-trait designs have been run with both verified software versions.

Fixed effects can apply to all traits or a named subset through `fixed[].traits`. Categorical effects and continuous covariates are supported. Additive genetic effects apply to every trait. Each additional random term has a unique name, a role (`litter`, `pe` or `iid`) and its own applicable traits. Covariances may be estimated between traits within a term; different named random terms are independent.

Litter effects need a litter ID, or the dam ID together with offspring birth date. Permanent environment is grouped by animal. Repeated records need an animal-and-event record key. For mixed designs, each trait needs an explicit event filter: keep a growth measurement on its original row and each farrowing on its own row, rather than copying growth phenotypes into farrowing records. Residual covariances between traits with no paired observations need explicit zero constraints; the workflow does not choose these assumptions automatically.

The current scope is Gaussian pedigree models. Count traits are checked for non-negative integers and fitted with a Gaussian model. Genomic, maternal genetic, heterogeneous-residual, threshold, random-regression and unknown-parent-group models are outside the current scope. Covariance between different named random terms is not supported.

### Joint growth and reproduction example

The [multitrait example](skills/pedigree-genetic-evaluation/references/example/multitrait-project.json) combines age at 120 kg and backfat measured once with repeated total-born records. Sex applies to the two growth traits, while parity and farrowing season apply to total born. The `birth_litter` term applies to growth and `sow_pe` to total born. Additive genetic covariance links all three traits; growth–total-born residual covariances are explicitly fixed at zero because the observations come from different events.

Copy the specification before editing it. From the repository root in PowerShell:

```powershell
Copy-Item skills/pedigree-genetic-evaluation/references/example/multitrait-project.json skills/pedigree-genetic-evaluation/references/example/multitrait-local.json
```

Check that copy without starting a solver:

```sh
Rscript --vanilla skills/pedigree-genetic-evaluation/scripts/run_eval.R --spec skills/pedigree-genetic-evaluation/references/example/multitrait-local.json --engine blupf90 --prepare-only
```

After software setup, remove `--prepare-only` for a full run. Use `--engine dmu` or `--engine compare` for DMU or both programs. Keep data paths relative to the copied JSON, or change them to absolute paths if moving it to another directory. The original [two-analysis example](skills/pedigree-genetic-evaluation/references/example/project.json) remains available for separate growth and reproduction analyses.

## Inputs and results

The project JSON specifies files, columns, traits, units, model and QC. Excel, CSV, TSV and delimited text inputs are supported. The [bundled example](skills/pedigree-genetic-evaluation/references/example/project.json) contains simulated data only; see its [provenance](skills/pedigree-genetic-evaluation/references/example/provenance.md).

Default QC applies mean ±3 SD once, separately for each trait. A missing or rejected observation does not remove another valid trait on the same animal. Parent information can be supplemented from other files; conflicting relationships and duplicate events are reported for review.

| Output | Contents |
| --- | --- |
| `delivery/evaluation_summary.xlsx` | QC, descriptive statistics and CV, variance components, heritability, repeatability or genetic correlation where applicable, standard errors and annual summaries |
| `delivery/analysis_*/EBV.csv` | Original animal IDs, EBV, SEP and reliability where available |
| `delivery/analysis_*/*.png` | QC histograms and annual trend plots |
| Run directory | Input snapshots, reproduction configuration, code, parameter files, ID maps, solver logs, versions and run manifest |

Variance components are reported separately for each named random term. A trait's phenotypic variance includes all terms acting on it; repeatability is reported only when a permanent-environment term applies, with additive and permanent-environment variance in its numerator. Litter and other independent grouping effects are not included in that numerator.

Standard errors for derived genetic parameters use the full sampling covariance matrix and the Delta method. Annual summaries use measurement or farrowing year by default; each animal contributes one EBV per trait and year. In a mixed design, only the event rows assigned to that trait are used, so a sow's farrowing years do not become years of growth evaluation. See [result definitions](skills/pedigree-genetic-evaluation/references/results.md) for interpretation and numerical checks.

Input paths and a relative `--output` are resolved from the **project JSON directory**. The bundled example therefore writes to `skills/pedigree-genetic-evaluation/references/example/example_runs/`, with a separate subdirectory for each run.

Add `--resume <run-directory>` to continue an existing run. Matching data, model and program fingerprints allow verified results to be reused. Changing only annual summaries or plots does not require another variance-component fit.

## Documentation

- [Project configuration](skills/pedigree-genetic-evaluation/references/project-spec.md)
- [BLUPF90 adapter](skills/pedigree-genetic-evaluation/references/blupf90.md)
- [DMU adapter](skills/pedigree-genetic-evaluation/references/dmu.md)
- [Results and standard errors](skills/pedigree-genetic-evaluation/references/results.md)
- [Skill instructions](skills/pedigree-genetic-evaluation/SKILL.md)

Every guide has an English and Simplified Chinese version, linked at the top of the page.

## Tests

Run the offline data-processing and statistical tests without either solver:

```sh
Rscript --vanilla tests/run_tests.R
```

Test fixtures and the bundled examples are synthetic. The offline suite covers data preparation and statistical calculations, including trait-specific effects, mixed single/repeated records and multiple random terms. Native REML and BLUP checks cover the original models, three- and four-trait designs, partial trait missingness, inbreeding and cross-program comparisons with the verified versions listed above.

To run the native regression checks, configure both programs in the local profile, then run:

```sh
Rscript --vanilla tests/run_native_tests.R --profile skills/pedigree-genetic-evaluation/host_profile.json
```

This runs both bundled project specifications in comparison mode and checks that resuming the completed runs does not restart the solvers. Add `--output <absolute-directory>` to choose where the test results are saved.
