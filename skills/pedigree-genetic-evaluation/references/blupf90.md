[English](blupf90.md) | [简体中文](blupf90.zh-CN.md)

# BLUPF90 backend: versions and file formats

The adapter supports RENUMF90 **1.169** and BLUPF90+ **2.76**. Executables are supplied through the host software configuration and are not distributed with the skill. The version is read from actual stdout. Other versions are rejected until their output formats have been tested; a newer release is not assumed to be compatible.

## Estimation and prediction

`f90_fit(ctx, directory)` runs renumbering and `OPTION method VCE` in a separate attempt directory to estimate variance components. `f90_blup(ctx, fit, directory)` writes the final fitted components to another parameter file and predicts EBV and prediction standard errors using `OPTION method BLUP`, `OPTION solv_method FSPAK`, and `OPTION sol se`. The VCE `solutions` file is not used as the final EBV/SEP source.

Each attempt preserves data, pedigree, both parameter files, stdout, stderr, the model-column map, and status JSON. Directories containing prior program files are not overwritten. Reuse is governed by input, model and software fingerprints. A `FINISHED` message alone does not establish that the result belongs to the current input.

## Data columns and model effects

The exported `phenotypes.txt` is whitespace-delimited without a header. Columns are ordered as traits `y1 ... yN`, all fixed-effect columns, the public integer animal `id`, and the additional random-group columns. Shortcuts retain the optional `shared` column; general models use separate grouping columns for named terms. `f90_model_map.json` records the actual column and effect numbers. Trait order always comes from the common model specification, not from the magnitude of the results.

Fixed classifications are exported as `cross alpha`; numeric covariates as `cov`, without an additional alpha/numer qualifier. The explicit mean is retained as a fixed classification. Each effect's position vector has a data-column number for applicable traits and 0 for inactive traits. Fixed effects can therefore differ among traits; a zero position means the effect is absent, not that the observation is zero.

The additive effect number is the number of fixed effects plus 1; additional terms follow in their declared order. Effect 3 must therefore not be hard-coded as additive, or effect 4 as litter. `RANDOM_GROUP` refers to an effect row, whereas observation and effect positions refer to data columns. Native effect numbers are mapped to unique component names such as `birth_litter` or `sow_pe`.

Common litter, permanent environment and generic grouped `iid` terms are exported as separate `diagonal` random effects. The model specification records each component's name, semantic role and active trait subset; the program's `Genetic variance(s)` heading does not determine scientific interpretation. A PE code must equal the public animal code on applicable observed rows, linking all litters from one sow to the same PE level. Different named random terms are independent; covariances are estimated within a term's active trait block, not between separate terms.

There is no fixed two-trait limit. Single and repeated traits can be combined with explicit event membership, trait-specific effects and justified residual-zero constraints. The [general-model specification](project-spec.md) explains these checks; the [mixed-record example](example/multitrait-project.json) gives a complete three-trait input. More traits do not remove identifiability or runtime limits.

The public pedigree is the validated closure of target animals and their ancestors. Its three exported columns are animal, sire and dam, with unknown parents coded 0. `PED_DEPTH 0` retains the complete selected closure, including target animals without their own records, and uses a complete ancestor search. Trimming that pedigree inside the backend would change the common basis of the two software models and is not permitted.

## Relationship-matrix precision

In RENUMF90 1.169, column 4 of `renaddXX.ped` stores the inbreeding-adjusted Mendelian sampling inverse as an integer, at approximately 0.001 resolution. In complex inbred pedigrees, this rounding can accumulate in A-inverse and affect EBV/PEV comparisons even when the same variance components are fixed. Relationship-matrix precision must therefore be checked rather than merely widening result tolerances.

With inbreeding enabled, the adapter reads `renf90.inb`, verifies F and both animal numberings, and constructs a higher-precision A-inverse using Henderson updates. Each animal contributes `Dinv × v × v'`, where v is 1 for the animal and −0.5 for each known parent, and `Dinv=4/(4−number of known parents−sire F−dam F)`. R's **Matrix** package combines sparse entries. The upper triangle is written as i/j/value with 17 significant digits.

VCE and final BLUP use the same file through the supported `RANDOM_TYPE user_file` interface. Native renumbered pedigrees are retained alongside `pedigree_Ainv_high_precision.txt` and `pedigree_Ainv_precision.json`. The JSON records numbering, construction method, native rounding differences and precision sources. This path depends only on F90 and R; it does not call DMU.

The source F values still have six decimal places in the program output; they are **not mathematically exact inbreeding coefficients**. This matrix reduces the integer-encoding error while retaining finite input precision. Fitting and BLUP must use the same precision. Cached estimates from the older native-integer path are explicitly rejected and require refitting.

With inbreeding disabled, the native F90 non-inbreeding path is retained and F is returned as NA with the setting recorded. Inbreeding is not silently disabled, parents are not changed, and DMU is not introduced to remove software differences. The actual `user_file` type in the software input and logs must remain visible in the delivery record. See the [official F90 user_file documentation](https://nce.ads.uga.edu/wiki/doku.php?id=user_defined_files_for_covariances_of_random_effects).

## Missing values and covariance constraints

Internal data use NA; F90 files use **-999**, with explicit `OPTION missing -999`. An actual observation equal to -999 is rejected to avoid collision with the missing marker. A record is retained when at least one trait is valid; records missing every trait cannot enter the backend. IDs must be present, while fixed effects and random groups are required on their applicable observed rows. Neutral values on inactive rows are input placeholders, not new observations.

BLUPF90+ VCE treats an initial 0 as a fixed 0. Unconstrained covariance starts must therefore be nonzero, and fixed zeros must be declared in `zero_covariances` in the project specification. Each component is a full global N×N matrix. Rows and columns outside its active traits are structural zeros; its active principal submatrix must be symmetric and positive definite. Matrices are written with 17 significant digits. Fixed-zero covariances are supported, not implicit removal of an active random effect through a zero variance. Remove such an effect explicitly from the model.

Trait pairs with no real co-observation require an explicit residual-zero constraint before either backend runs. For identity-related terms, pairs that never share a group likewise require a justified zero covariance. User and structural zeros are excluded from estimated-parameter metadata; the component name and global trait indices determine every free parameter.

RENUM may print starts in a shorter format. Before VCE, the adapter rewrites the original matrices at full precision so that a small nonzero covariance is not rounded to zero. Final fitted matrices are parsed from actual program output, not the original parameter file or a manually entered table. The 2.76 text output generally provides about five significant digits. All available digits are retained without further two-decimal rounding; precision absent from native output cannot be recovered.

## IDs, inbreeding and prediction errors

RENUM renumbers the public integer IDs again. The verified `renaddXX.ped` format has 10 columns: column 1 is the F90 internal code, and column 10 is the input public ID. The adapter checks column count, a one-to-one map and complete public-pedigree coverage, then recovers original string IDs through the public pedigree. Output `id` is the public code; F90 internal numbers cannot be used to merge separate models.

In long-format EBV output, `has_record` means that an animal has at least one nonmissing own observation for that trait. An animal with age records only has TRUE for age and FALSE for backfat. It may still receive a backfat EBV through the multitrait model and relatives. A male without TNB observations can similarly receive a TNB EBV from pedigree and correlated-trait information. Overall model membership is recorded separately through public data animal codes.

Inbreeding adjustment is required by default. In the verified version, `renf90.inb` has three columns: public ID, F and F90 internal code. Both numberings must be checked before reading F. An anonymous full-sib mating case verifies F=0.25. If verified F is unavailable, the adapter returns NA and a reason; unknown inbreeding is never treated as F=0. With adjustment disabled, F is returned as NA and the setting is recorded.

Fixed-component BLUP `solutions` must have five columns: trait, effect, internal level, solution and SE. Only the current additive effect is extracted, with exactly one result required for every public animal-by-trait combination. Under `OPTION sol se`, this random effect's SE is the **prediction standard error SEP**, and PEV=SEP². Fixed-effect SE, variance-component SE and SEP are not interchangeable.

BLUP also requests `OPTION residual`. In 2.76 `yhat_residual`, the first n trait columns are fitted values, followed by n residual columns. The adapter verifies renumbered animal/record order, missing masks, and fitted plus residual equals the original observation. It then attaches `attr(ebv, "fitted")` with columns record_id/trait/predicted for fixed-component comparisons. Native fitted zeros at missing observations are converted to NA rather than treated as predictions.

## Sampling covariance and convergence

The complete upper triangle of `AI_inv.txt` is preferred. Its seven columns in 2.76 are sampling-matrix row, column, the row parameter's effect number, trait1, trait2, group dimension and sampling covariance. Diagonal records provide each parameter's role/trait metadata. The residual effect code is the total number of effects plus 1, with group dimension 0; random groups have dimension equal to the trait count.

From these fields, the adapter builds `param_meta` and checks every estimated parameter, upper-triangle coverage, symmetry and positive definiteness. Its `role` field identifies the unique component key; the model's random-term registry supplies its semantic role. Trait indices use the global trait order, including for nonconsecutive subsets. Fixed-zero and inactive parameters are absent from the sampling matrix. If the file is unavailable, the complete printed inverse-AI matrix from the same verified version may be used with a determined parameter order. Invalid format or an unreliable matrix yields `param_cov=NULL` with a reason. Squared component SEs are never assembled into a diagonal substitute.

`native_se` retains native SE matrices under the same keys as the fitted components: additive, each declared random term, and residual. Shortcuts retain `shared`. These can support individual component SEs but cannot replace the full sampling covariance needed for heritability, repeatability or genetic-correlation SEs. The shared results module uses the full parameter covariance and the delta method.

Exit code, completion marker, final convergence criteria, final matrices and AI information are checked separately. Boundaries, singularity, bending and NaN are recorded as diagnostics. Program completion and model validation are distinct states. On Windows, a missing parameter file can still produce exit code 0, so generated files and actual completion information must also be checked.

## Cross-checks

The repository's synthetic examples retain the original bivariate litter and single-trait repeated-record projects and add a joint growth/reproduction model. Native checks also cover three/four traits, nonconsecutive active subsets and complete sampling-covariance mapping. With actual executable paths configured, examples can check the full workflow. Public offline tests cover preprocessing and result calculations; they do not replace solver runs.

On Windows, `processx` sets environment variables without inserting them into parameter files or program arguments. Each process preserves stdout, stderr and exit status.

Cross-checks account for native file precision: DMU SOL has approximately six significant digits, F90 solutions eight fixed decimal places, and native F90 `yhat_residual` approximately eight significant digits; internal calculations can also have differences at single-precision scale. PEV tolerance propagates SEP rounding bounds. The numerical allowance for fitted values uses the prediction scale, not an unusually small phenotype as a substitute. These rules are defined before execution, not relaxed after seeing a result.

Primary references: [RENUMF90 input and renumbered files](https://nce.ads.uga.edu/wiki/doku.php?id=readme.renumf90), [BLUPF90+ modes and sol se](https://nce.ads.uga.edu/wiki/doku.php?id=readme.blupf90plus), and [zero-covariance constraints in the F90 tutorial](https://nce.ads.uga.edu/wiki/lib/exe/fetch.php?media=tutorial_blupf90.pdf). Format details were also verified with actual anonymous runs of versions 1.169/2.76.
