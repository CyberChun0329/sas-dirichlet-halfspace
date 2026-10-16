# Paper-to-evidence map

All paths below are repository-relative. Let `V` denote
`runs/correctness_01/work/results`, `X` denote `runs/correctness_extra_01`,
and `A` denote `runs/paper_minimum_aggregate_02`.

| Paper object | Saved source | Fortran producer |
|---|---|---|
| Figure 1; Case I validation; Table S1 | V/case1_profiles.csv, case1_order_check.csv, case1_summary.txt; X/case1_z_ge_half.csv | validation_driver; stage_c_reduce |
| Table 2: 50-depth accumulated deficit performance | A/cpu_deficit_comparisons.csv; raw task deficit_original_0001_N1 | targeted_driver; replay_summary |
| Figure 2; Table 3; Table S4 | V/case2_profiles.csv, case2_convergence.csv; X/gaussian_coeff_*; runs/original_gaussian_controls_01 | validation_driver; gaussian_coeff_check; gaussian_original_controls |
| Table 4; Figure 3 | V/case3_table.csv, case3_grid_1…8.csv; X/case3_six_reductions.csv, case3_z4_ratios.csv | validation_driver; stage_c_reduce |
| Table 5; Figure 4; main performance scope | A/cpu_main100_comparisons.csv; paper_0001…0100_N1 raw batches | targeted_driver; replay_summary |
| Figure 5 | A/cpu_budget3_budgets.csv; B1 paper_0076 and paper_0080 | targeted_driver; replay_summary |
| Tables S5–S6; supplementary periodic response and history bounds | V/periodic_grid_1…8.csv, periodic_refinement_union.csv, periodic_summary.txt; X/periodic_* | periodic_driver; stage_c_reduce |
| Table S2 | V/table_s2.csv | validation_driver |
| Table S3; validation/precision counts | V/stability.csv, core_summary.txt, precision.csv, precision_summary.txt and companion regressions | validation_driver; precision_driver; summarize_validation |
| Table S7, preparation and N=1/10/100 | A/cpu_repeat10_comparisons.csv, cpu_repeat100_comparisons.csv, preparation.csv and matching N1 main rows | targeted_driver; replay_summary |
| Table S8; stronger spatial routes | A/cpu_shared7_comparisons.csv | targeted_driver; replay_summary |
| Figure S1; deep/long-time reversal | A/cpu_deep_Q_comparisons.csv, cpu_deep_Q_budgets.csv | targeted_driver; replay_summary |
| Deep uniform diffusion-scale comparisons | A/cpu_uniform_N1_comparisons.csv; repeat batches | targeted_driver; replay_summary |

`provenance/figure_sources.json` records the plotting rows, fields and hashes
using repository-relative paths. The source CSVs are under `runs/`, and the
distributed figures are in `figures/`. Regenerated images go to `build/`.

## Scope and interpretation

The main comparator is the preprocessed, analytically time-integrated
**two-dimensional** spatial kernel. Seven selected cases also measure SHARED_1D;
ALL_SPACE is only the union of actually measured non-SAS keys in those cases,
not all possible spatial algorithms or a 100-case one-dimensional comparison.

ACTUAL_ERROR pairs match same-order **observed** error; they are not the same as
TARGET pairs under a loose error cap. QUERY reuses prepared contexts. TOTAL
prepares and evaluates in its own timed callback; it is not QUERY plus an
independently measured PREP estimate. CPU and wall observations stay separate.
Original formal batches ran with up to eight concurrent processes and one-thread
environment settings. Timing results describe that execution environment.

Budget traces describe only the explicitly measured finite candidate lists.
`runs/paper_minimum_inputs_01/shortlist_diagnostics.csv` and the remaining-input
diagnostics retain omitted intermediate budget nodes. Reporting-floor segments,
unavailable intervals, finite endpoints, conservative directions and crossings
remain relevant; selected large gains are not continuous optimal-frontier bounds.

The independent validation kernel is separate from SAS. Performance references
are independent of SAS but share implementation with some spatial routes; their
refinement floors are empirical, not rigorous true-error certificates. Periodic
extrema are numerically refined, not interval-certified. The periodic validation
covers 96 depths; additional refinement covers the 12 table depths and two depths
at L=0.2, Q, d=0.1. The maxima over that subset refer to all 14 depths.

The periodic CSVs retain the recorded temperature and refinement fields in their
original row order. `provenance/source_records.json` records their source hashes
and column selections. It also identifies the saved byte-equality checks for the
16 distributed CPU/wall summary reports. The Fortran `replay_summary` source and
its saved inputs are included; those stored checks are distinct from fresh
numerical execution on another machine.
