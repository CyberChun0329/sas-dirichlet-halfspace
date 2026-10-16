# Fortran replay summary

Inputs: runs/paper_minimum_aggregate_02/cpu_repeat100_comparisons.csv and runs/paper_minimum_aggregate_02/cpu_repeat100_budgets.csv
Comparison CSV rows: 20
Distinct comparison case IDs: 5
Absolute-only evidence rows: 0
Budget CSV rows: 0
Distinct budget case IDs: 0

Counts describe exact saved rows and distinct case IDs within each group. Cases recur across families, modes, caps and budget directions; group counts must not be pooled as independent cases.
No population confidence intervals are inferred. Speed ranges describe recorded cost ratios. A cost winner can remain at the reporting floor; it is not a resolved precision claim.

## DIRECT_2D QUERY ACTUAL_ERROR cap=  1.0000000E-06
rows=5 distinct_case_ids=5
unmatched=0 unreached=0 floor_pairs=0 resolved_pairs=5 one_floor_pairs=0
SAS_cost_wins=3 SPACE_cost_wins=2 resolved_SAS_wins=3 resolved_SPACE_wins=2 ties_or_vetoed=0
separated_winners=5 rows_with_dominance_veto=0
SPACE_dominates_SAS_rows=2 SAS_dominates_SPACE_rows=3
SPACE/SAS cost ratio (selected pairs):   6.0434622E-01 to   4.6374023E+01
SAS winner speed ratio:   3.7415737E+01 to   4.6374023E+01
SPACE winner speed ratio:   1.5025359E+00 to   1.6546806E+00
status ACTUAL_ERROR_MATCH: 5

## DIRECT_2D QUERY TARGET cap=  1.0000000E-06
rows=5 distinct_case_ids=5
unmatched=0 unreached=0 floor_pairs=0 resolved_pairs=5 one_floor_pairs=0
SAS_cost_wins=3 SPACE_cost_wins=2 resolved_SAS_wins=3 resolved_SPACE_wins=2 ties_or_vetoed=0
separated_winners=5 rows_with_dominance_veto=0
SPACE_dominates_SAS_rows=0 SAS_dominates_SPACE_rows=0
SPACE/SAS cost ratio (selected pairs):   6.0434622E-01 to   4.6374023E+01
SAS winner speed ratio:   3.7415737E+01 to   4.6374023E+01
SPACE winner speed ratio:   1.5025359E+00 to   1.6546806E+00
status TARGET_BOTH_REACHED: 5

## DIRECT_2D TOTAL ACTUAL_ERROR cap=  1.0000000E-06
rows=5 distinct_case_ids=5
unmatched=0 unreached=0 floor_pairs=0 resolved_pairs=5 one_floor_pairs=0
SAS_cost_wins=3 SPACE_cost_wins=2 resolved_SAS_wins=3 resolved_SPACE_wins=2 ties_or_vetoed=0
separated_winners=5 rows_with_dominance_veto=0
SPACE_dominates_SAS_rows=2 SAS_dominates_SPACE_rows=3
SPACE/SAS cost ratio (selected pairs):   6.0146843E-01 to   1.0187174E+02
SAS winner speed ratio:   6.7215600E+01 to   1.0187174E+02
SPACE winner speed ratio:   1.5003342E+00 to   1.6625977E+00
status ACTUAL_ERROR_MATCH: 5

## DIRECT_2D TOTAL TARGET cap=  1.0000000E-06
rows=5 distinct_case_ids=5
unmatched=0 unreached=0 floor_pairs=0 resolved_pairs=5 one_floor_pairs=0
SAS_cost_wins=3 SPACE_cost_wins=2 resolved_SAS_wins=3 resolved_SPACE_wins=2 ties_or_vetoed=0
separated_winners=5 rows_with_dominance_veto=0
SPACE_dominates_SAS_rows=0 SAS_dominates_SPACE_rows=0
SPACE/SAS cost ratio (selected pairs):   6.0146843E-01 to   1.0187174E+02
SAS winner speed ratio:   6.7215600E+01 to   1.0187174E+02
SPACE winner speed ratio:   1.5003342E+00 to   1.6625977E+00
status TARGET_BOTH_REACHED: 5

Budget gains count only intervals with end > start: 10-to-100 means >=10 and <100; 100 means >=100. Finite terminal endpoints and zero-width rows are retained separately.
