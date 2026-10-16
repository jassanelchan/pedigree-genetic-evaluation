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

The repository provides R scripts and is not currently an installable R package.

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
| Animal | 1 or 2 | Additive genetic, residual |
| Common litter | 1 or 2 | Additive genetic, common litter, residual |
| Repeatability | 1 | Additive genetic, permanent environment, residual |

Both traits use the same fixed-effect structure, with categorical effects and continuous covariates supported. Litter models need a litter ID or the dam ID together with offspring birth date. Repeated records need an animal-and-event record key.

The current scope is Gaussian pedigree models. Count traits are checked for non-negative integers and fitted with a Gaussian model. Genomic, three-trait, threshold, random-regression, maternal genetic and unknown-parent-group models are outside the current scope.

## Inputs and results

The project JSON specifies files, columns, traits, units, model and QC. Excel, CSV, TSV and delimited text inputs are supported. The [bundled example](skills/pedigree-genetic-evaluation/references/example/project.json) contains simulated data only; see its [provenance](skills/pedigree-genetic-evaluation/references/example/provenance.md).

Default QC applies mean ±3 SD once, separately for each trait. A missing or rejected observation does not remove another valid trait on the same animal. Parent information can be supplemented from other files; conflicting relationships and duplicate events are reported for review.

| Output | Contents |
| --- | --- |
| `delivery/evaluation_summary.xlsx` | QC, descriptive statistics and CV, variance components, heritability, repeatability or genetic correlation where applicable, standard errors and annual summaries |
| `delivery/analysis_*/EBV.csv` | Original animal IDs, EBV, SEP and reliability where available |
| `delivery/analysis_*/*.png` | QC histograms and annual trend plots |
| Run directory | Input snapshots, reproduction configuration, code, parameter files, ID maps, solver logs, versions and run manifest |

Standard errors for derived genetic parameters use the full sampling covariance matrix and the Delta method. Annual summaries use measurement or farrowing year by default; each animal contributes one EBV per year. See [result definitions](skills/pedigree-genetic-evaluation/references/results.md) for interpretation and numerical checks.

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

Test fixtures and the bundled example are synthetic. The adapters have also been checked with native REML and BLUP runs covering partial trait missingness, inbreeding, common litter, permanent environment and cross-program comparisons.
