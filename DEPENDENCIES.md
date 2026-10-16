# Build dependencies and validation boundary

Required: GNU Fortran with `real128` support, GNU Make, MPFR and GMP development
libraries, and a BLAS/LAPACK implementation. The checked platform is Apple Silicon
macOS, GNU Fortran 15.2.0, using Apple Accelerate and Homebrew MPFR/GMP. Recorded
compiler and library identities are in `provenance/source_records.json`.

The Makefile defaults to `/opt/homebrew/opt/mpfr` and `/opt/homebrew/opt/gmp`.
Override `MPFR_PREFIX`, `GMP_PREFIX`, `MPFR_LIBS` and `ACCELERATE` when required.
For example, a non-macOS BLAS/LAPACK installation would need an appropriate
`ACCELERATE='-llapack -lblas'`; that port has not been tested here.

```sh
make -j2 release          # 9 executables; compile only
make -j2 debug            # checked diagnostics build; compile only
make -j2 build-reference-generators  # opt-in polynomial reference generator
make -j2 build-tests      # optional; compile tests, do not run them
```

There is no automatic dependency installation. All products stay below `build/`,
which is Git-ignored. Module output directories are separate for parallel builds.
The default release flags retain `-O3 -march=native`, not `-ffast-math`.

Native reference streams are compiler/platform representations with strict
input/release checks, not a portable interchange format. Different flags or
compilers can cause identity rejection; do not replace checks with tolerances.
Debug builds previously detected last-bit coefficient differences against
release-created polynomial references. Readable CSV evidence remains available;
fresh native reference generation on another platform requires its own validation.

Node 22 or newer is used for file verification and the explicit runner. Optional
figure rendering needs Python 3 and Matplotlib (`tools/render_figures.py`).
That script reads saved CSV values and writes `build/rendered_figures/`; it does
not solve the heat equation or recompute errors, frontiers or statistics.

The recorded build checks on the platform above covered release/debug
executables, optional reference generators and test executables. Numerical
execution and timing records are indexed separately in `PAPER_RESULTS.md`.
