[English](dmu.md) | [简体中文](dmu.zh-CN.md)

# DMU pedigree backend

## Supported runtime

The implemented adapter is for **DMU version 6, release 5.6**, using the
separate `dmu1`, `dmuai`, and `dmu4` executables. The release printed by the
actual `dmu1` process is checked. A different release or the newer single-DMU
layout is not silently treated as compatible.

Set `software.dmu_dir` to the verified `bin` directory (or its parent).
Executable binaries belong to the host profile, not inside this skill.
The adapter invokes the modules directly and preserves all staged input,
intermediate files, stdout, stderr, and process metadata. It does not use the
distributed batch wrappers, which delete intermediate files.

Official primary sources:

- [R5.6 user guide](https://dmu.ghpc.au.dk/dmu/DMU/Doc/Current/dmuv6_guide.5.6.pdf)
- [Windows packages](https://dmu.ghpc.au.dk/dmu/DMU/Windows/Current/)
- [Linux packages](https://dmu.ghpc.au.dk/dmu/DMU/Linux/Current/)
- [Official parameter-SE explanation](https://dmu.ghpc.au.dk/dmu/DMU/Doc/Current/STDERR.pdf)

The official newer HTML guide describes another executable/launcher layout
and conditional research-license output limits. Do not infer that those
limits apply to every R5.6 installation. Verify the installed version,
authorization, and actual target-animal output coverage. Research use should
acknowledge DMU. Commercial routine evaluation requires agreement with QGG.

## Runtime configuration

Obtain the compatible release from the official distribution and set
`software.dmu_dir` to its `bin` directory. Keep the release's shared runtime
libraries with the executables. The adapter runs the modules directly rather
than the batch wrappers, preserving intermediate files and process logs.

The runtime is supplied by the installation environment and kept outside the
skill directory. Sourcing the R modules does not download, install, or execute
third-party programs.

## Models and input contract

Supported models are Gaussian pedigree animal models with one or more traits,
trait-specific fixed effects and multiple independent identity-related random
terms. The animal, litter and repeatability shortcuts remain available;
`model:"general"` combines named litter, animal permanent-environment and `iid`
terms on chosen trait subsets. The explicit mean and numeric covariates are
retained. See the [model specification](project-spec.md) for record filters,
grouping and identifiability requirements. There is no fixed two-trait limit,
but the data and runtime must support the requested model.

Genomic models, UPG/metafounders, maternal genetic effects, heterogeneous
residuals, threshold models and random regression are not implemented.

The shared context is translated to DMU data with INTEGER*4 classification
columns first and REAL*4 trait/covariate columns afterward. These two column
lists have separate 1-based indices in the model directives. Public integer
animal codes are shared with the other backend; original string IDs are saved
in the common mapping. Class code zero is not a valid observed fixed level.

Each trait's model lines list only its active fixed classes, covariates and
random factors. Values on inactive rows are neutral placeholders, not imputed
phenotypes. Event filters distinguish single growth records from repeated
farrowing records; copying growth values onto farrowing rows is invalid.

Missing real values use an engine-specific sentinel strictly below an
engine-specific threshold, selected below all finite staged values. The actual
DMU1 per-trait effective observation counts must match the original NA masks.
No extra QC is performed by this backend and no missing trait is fabricated.

`pedigree.txt` contains public animal code, sire, dam, and validated
parent-before-offspring sort. Unknown parents are zero. `PED 1` is used when
inbreeding is enabled and `PED 2` when it is disabled. All provided pedigree
codes must be unique and the ancestor closure complete.
With `PED 2`, the non-inbreeding relationship model is a modeling choice,
not evidence that every animal's actual F is zero. The adapter returns F as
NA, records that it is unverified, and does not allow a verified reliability
calculation from this mode.

Random factor 1 is always additive; named identity-related terms follow in
declaration order, then the residual. A PE is one-to-one with animal identity
and connects its repeated records; it is not a birth litter. Each factor has
its own active trait list. Its native covariance indices are local to that
list, so native index 2 may refer to global trait 3 when the subset is [1,3].
The adapter maps factor numbers and local indices back to unique component
names and global trait indices. Each component returned to the common workflow
is a full N×N matrix with zeros outside its active principal block. Different
named terms are independent; no between-term covariance is fitted.

Free covariance entries and structural-zero constraints are separate. In this
release a fixed-zero restriction is written with a comma in the pair, e.g.
`$VAR_REST COV 2 E 2,1`. DMU1 can print semantic errors while returning exit zero;
the adapter rejects its ERROR/ERRORS messages before starting the next module.
For residual pairs with no real co-observation, the common specification must
declare a zero constraint; the adapter writes DMU's non-existing-pair list.
It does not invent this assumption. Identity-related trait pairs that never
share a group likewise need an explicit justified zero. Constrained values
are checked against PAROUT where present, and missing native entries are
accepted only for declared fixed zeros.

## Estimation, prediction, and bounded execution

REML uses task 1, method 1 (AI with EM recovery). The orchestrator may request
one final bounded recovery attempt with `stable_reml=TRUE`; that selects task 1,
method 4 (AI with step halving). This does not change the scientific model.
Independent attempt directories retain both starts and failures. The total
wall-clock budget for DMU1 and its subsequent module is `ctx$timeout`; the
common process runner terminates the process tree at the deadline. Thread
counts use `MKL_NUM_THREADS` and `OMP_NUM_THREADS`.

Numerical convergence is accepted only when native output reports it. R5.6 may
accept a Hessian convergence criterion even when its update or gradient norms
exceed their separate thresholds; this is recorded explicitly. Successful
process completion alone is insufficient.

BLUP uses task 11, DMU4 direct FSPAK method 9 with the final components. The SOL
descriptor must show seven integer and two real fields. Only additive-factor
random solutions are exported as EBV. The original public animal code is the
fifth integer field; the two real fields are EBV and SEP for this verified
solver. DMU5's second real field is a previous-iteration solution and is never
interpreted as SEP. Exact animal-by-trait output coverage is mandatory.

The verified Windows R5.6 INBREED file prints `1+F` in column 8, despite the
manual's abbreviated field description. Founder values verify this scale;
the adapter subtracts one only when founder scaling supports it. A zero-based
F variant can also be recognized from founders. Other layouts are rejected.
The public INBREED text is rounded (approximately 0.001 in the verified file),
so returned F can have about 0.0005 rounding error. This does not change the
native relationship matrix used for solving.

Requested native predictions are parsed from R5.6 RESIDUAL: record number,
observed, predicted, residual, hat diagonal, and observed corrected for fixed
effects (one block per trait). Original missing observations remain NA in
the returned fitted table. For present observations, predicted plus residual
must reproduce the staged observation within REAL*4/print precision.

## Components, sampling covariance, and reporting

PAROUT records identify factor, row, column, and final component. PAROUT_STD
contains NPAR, parameter estimates with SE, and the complete lower-triangular
sampling correlation/covariance matrix. The R5.6 order is checked against
PAROUT. Restricted-zero parameters still occur as zero rows/columns in native
PAROUT_STD and are removed from the common *estimated-parameter* covariance.
Fixed-zero restrictions remain explicit in the model manifest. Native local
indices are mapped through each factor's active trait list before sampling
parameters are matched. The common `param_meta.role` field is a unique
component key; semantic `litter`, `pe` or `iid` roles come from the model registry.
Inactive structural zeros do not enter the estimated-parameter covariance.

Validation checks the full parameter-pair coverage, native SE squared against
sampling diagonals, correlation/covariance consistency, and positive
semidefiniteness within print precision. Missing or unreliable sampling
covariance is reported as unavailable; it is never replaced with a diagonal
matrix or SOL_STD, which concerns covariance of effect solutions.

The shared [results module](results.md) calculates h2, repeatability and additive rg from
declared model components with analytic delta-method gradients. Native
intra-class correlations are a cross-check, not an implicit denominator choice.
For trait t, Vp sums its additive, residual and all applicable named random
variances; h2 is VA/Vp. Repeatability is reported only for traits with PE and
uses `(VA + sum of applicable PE variances)/Vp`. Litter and generic `iid` terms
contribute to the denominator, not the repeatability numerator. The familiar
shortcut formulas remain VA/(VA+VC+VE) for litter-model h2, VA/(VA+VPE+VE) for
repeatability-model h2, and (VA+VPE)/(VA+VPE+VE) for repeatability.

DMU reports a restricted-likelihood criterion plus a data/fixed-effects
constant. It is preserved as such and not directly compared with another
engine's raw -2LogL or AIC. DMU REAL*4 data and printed SOL precision also mean
cross-engine tests require documented absolute and relative tolerances, not
byte-identical output or merely high EBV correlation.

## Example and tests

The original synthetic example exercises common litter, repeated records,
partial bivariate observations, fixed classes and a continuous covariate.
The [mixed-record example](example/multitrait-project.json) jointly analyses
two growth traits and repeated TNB with different effects and trait subsets.
Native tests also cover three/four traits, nonconsecutive subsets, named
random blocks and complete sampling-covariance mapping.
After supplying the executable paths, `engine:"compare"` runs both adapters
and checks mapped components and fixed-component EBV/PEV and fitted values.
The repository's offline tests check preprocessing and statistical functions
without invoking either native solver.

Real farm inputs and program binaries are not included in the repository.
