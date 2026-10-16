# Rectangular Dirichlet half-space: paper code and evidence

Companion to **Semi-analytical transient heat conduction from a rectangular
Dirichlet patch in a three-dimensional half-space**, Chunxiao Zheng.

This repository contains the Fortran implementation and saved numerical evidence
used by the current manuscript and supplementary material. Scientific evaluation,
references, error calculations, selection and numerical summaries use Fortran.
The optional Python script only renders figures from saved scalar CSV fields;
Node scripts only verify files and launch explicit jobs.

## Start without running experiments

```sh
node tools/verify_package.mjs
make -j2
node tools/run_paper.mjs --help
```

`make` compiles only. The checked toolchain is Apple Silicon/macOS with GNU
Fortran, MPFR/GMP and Accelerate; see [DEPENDENCIES.md](DEPENDENCIES.md).

## What is included

- Steady/transient polynomial moments, prepared SAS and spatial routes, and
  the independently implemented validation reference.
- Cases I–III, Gaussian truncation/evaluation checks and stability tests.
- Supplementary periodic responses and history-tail checks, including the
  96-depth input grid and 14-depth refinement subset.
- The explicit 100-case low-degree performance collection, repeated queries,
  accumulated-deficit benchmark, deep uniform and deep-Q controls.
- Seven measured one-dimensional spatial comparisons and their ALL_SPACE union,
  with the recorded comparison states.
- Native timing windows, separate preparation/query/total observations,
  selected correctness values/references, finite-budget records and all six
  manuscript/supplement figures.

[PAPER_RESULTS.md](PAPER_RESULTS.md) maps the paper objects to files and explains
the comparison limits. [REPRODUCING.md](REPRODUCING.md) gives explicit reproduction
commands. The runner executes one explicitly named validation chunk or timing task.

## Source and evidence identities

`config/source_release.txt` identifies the numerical release recorded in the
saved results and native references.
`provenance/package_files.json` and `SHA256SUMS` identify the current distributed
files. The polynomial rows and values retain the union of the explicit task
case/key pairs. `provenance/source_records.json` records the original identities
and literal selections used for the distributed evidence. Publication-source
changes remove unused experiment branches while retaining the paper's numerical
formulas and native-reference format. Periodic data contain temperature swings,
history bounds and refinement checks on the specified 96-depth grid and 14-depth
subset. [PAPER_RESULTS.md](PAPER_RESULTS.md) identifies the evidence used by each
paper object.

The author's code retains the existing MIT license; third-party material is not
relicensed by it.
