# Fortran replay summary

Inputs: runs/paper_minimum_aggregate_02/wall_deep_Q_comparisons.csv and runs/paper_minimum_aggregate_02/wall_deep_Q_budgets.csv
Comparison CSV rows: 4
Distinct comparison case IDs: 1
Absolute-only evidence rows: 0
Budget CSV rows: 69
Distinct budget case IDs: 1

Counts describe exact saved rows and distinct case IDs within each group. Cases recur across families, modes, caps and budget directions; group counts must not be pooled as independent cases.
No population confidence intervals are inferred. Speed ranges describe recorded cost ratios. A cost winner can remain at the reporting floor; it is not a resolved precision claim.

## DIRECT_2D QUERY ACTUAL_ERROR cap=  1.0000000E-06
rows=1 distinct_case_ids=1
unmatched=0 unreached=0 floor_pairs=0 resolved_pairs=1 one_floor_pairs=0
SAS_cost_wins=0 SPACE_cost_wins=1 resolved_SAS_wins=0 resolved_SPACE_wins=1 ties_or_vetoed=0
separated_winners=1 rows_with_dominance_veto=0
SPACE_dominates_SAS_rows=0 SAS_dominates_SPACE_rows=0
SPACE/SAS cost ratio (selected pairs):   7.7300780E-01 to   7.7300780E-01
SAS winner speed ratio: unavailable
SPACE winner speed ratio:   1.2936480E+00 to   1.2936480E+00
status ACTUAL_ERROR_MATCH: 1

## DIRECT_2D QUERY TARGET cap=  1.0000000E-06
rows=1 distinct_case_ids=1
unmatched=0 unreached=0 floor_pairs=0 resolved_pairs=1 one_floor_pairs=0
SAS_cost_wins=0 SPACE_cost_wins=1 resolved_SAS_wins=0 resolved_SPACE_wins=1 ties_or_vetoed=0
separated_winners=1 rows_with_dominance_veto=0
SPACE_dominates_SAS_rows=0 SAS_dominates_SPACE_rows=0
SPACE/SAS cost ratio (selected pairs):   7.7300780E-01 to   7.7300780E-01
SAS winner speed ratio: unavailable
SPACE winner speed ratio:   1.2936480E+00 to   1.2936480E+00
status TARGET_BOTH_REACHED: 1

## DIRECT_2D TOTAL ACTUAL_ERROR cap=  1.0000000E-06
rows=1 distinct_case_ids=1
unmatched=0 unreached=0 floor_pairs=0 resolved_pairs=1 one_floor_pairs=0
SAS_cost_wins=0 SPACE_cost_wins=1 resolved_SAS_wins=0 resolved_SPACE_wins=1 ties_or_vetoed=0
separated_winners=1 rows_with_dominance_veto=0
SPACE_dominates_SAS_rows=0 SAS_dominates_SPACE_rows=0
SPACE/SAS cost ratio (selected pairs):   8.2515933E-01 to   8.2515933E-01
SAS winner speed ratio: unavailable
SPACE winner speed ratio:   1.2118872E+00 to   1.2118872E+00
status ACTUAL_ERROR_MATCH: 1

## DIRECT_2D TOTAL TARGET cap=  1.0000000E-06
rows=1 distinct_case_ids=1
unmatched=0 unreached=0 floor_pairs=0 resolved_pairs=1 one_floor_pairs=0
SAS_cost_wins=0 SPACE_cost_wins=1 resolved_SAS_wins=0 resolved_SPACE_wins=1 ties_or_vetoed=0
separated_winners=1 rows_with_dominance_veto=0
SPACE_dominates_SAS_rows=0 SAS_dominates_SPACE_rows=0
SPACE/SAS cost ratio (selected pairs):   8.2515933E-01 to   8.2515933E-01
SAS winner speed ratio: unavailable
SPACE winner speed ratio:   1.2118872E+00 to   1.2118872E+00
status TARGET_BOTH_REACHED: 1

Budget gains count only intervals with end > start: 10-to-100 means >=10 and <100; 100 means >=100. Finite terminal endpoints and zero-width rows are retained separately.

## DIRECT_2D QUERY nominal budget
rows=12 distinct_case_ids=1 positive_width_intervals=12 zero_width_rows=0 finite_endpoints=1
positive_width_SAS_10_to_100=1 SAS_100=1 SPACE_10_to_100=3 SPACE_100=1
recorded finite bank budget span:   3.9564787E-06 to   1.5770000E-03
classification UNAVAILABLE: 3
classification NEITHER_100: 6
classification SPACE_100: 1
classification SAS_100: 1
classification BOTH_RELATIVE_FLOOR: 1
all-width gain label NONE: 6
all-width gain label SPACE_10_to_100: 3
all-width gain label SPACE_100: 1
all-width gain label SAS_100: 1
all-width gain label SAS_10_to_100: 1

## DIRECT_2D QUERY SAS_conservative budget
rows=11 distinct_case_ids=1 positive_width_intervals=11 zero_width_rows=0 finite_endpoints=1
positive_width_SAS_10_to_100=1 SAS_100=1 SPACE_10_to_100=3 SPACE_100=2
recorded finite bank budget span:   3.9564787E-06 to   1.5770000E-03
classification UNAVAILABLE: 2
classification NEITHER_100: 5
classification SPACE_100: 2
classification SAS_100: 1
classification BOTH_RELATIVE_FLOOR: 1
all-width gain label NONE: 4
all-width gain label SPACE_10_to_100: 3
all-width gain label SPACE_100: 2
all-width gain label SAS_100: 1
all-width gain label SAS_10_to_100: 1

## DIRECT_2D QUERY SPACE_conservative budget
rows=12 distinct_case_ids=1 positive_width_intervals=12 zero_width_rows=0 finite_endpoints=1
positive_width_SAS_10_to_100=1 SAS_100=2 SPACE_10_to_100=3 SPACE_100=0
recorded finite bank budget span:   3.9564787E-06 to   1.5770000E-03
classification UNAVAILABLE: 3
classification NEITHER_100: 6
classification SAS_100: 2
classification BOTH_RELATIVE_FLOOR: 1
all-width gain label NONE: 6
all-width gain label SPACE_10_to_100: 3
all-width gain label SAS_100: 2
all-width gain label SAS_10_to_100: 1

## DIRECT_2D TOTAL nominal budget
rows=12 distinct_case_ids=1 positive_width_intervals=12 zero_width_rows=0 finite_endpoints=1
positive_width_SAS_10_to_100=1 SAS_100=1 SPACE_10_to_100=3 SPACE_100=1
recorded finite bank budget span:   4.5209040E-06 to   1.9446000E-03
classification UNAVAILABLE: 3
classification NEITHER_100: 6
classification SPACE_100: 1
classification SAS_100: 1
classification BOTH_RELATIVE_FLOOR: 1
all-width gain label NONE: 6
all-width gain label SPACE_10_to_100: 3
all-width gain label SPACE_100: 1
all-width gain label SAS_100: 1
all-width gain label SAS_10_to_100: 1

## DIRECT_2D TOTAL SAS_conservative budget
rows=10 distinct_case_ids=1 positive_width_intervals=10 zero_width_rows=0 finite_endpoints=1
positive_width_SAS_10_to_100=1 SAS_100=1 SPACE_10_to_100=2 SPACE_100=1
recorded finite bank budget span:   4.5209040E-06 to   1.9446000E-03
classification UNAVAILABLE: 2
classification NEITHER_100: 5
classification SPACE_100: 1
classification SAS_100: 1
classification BOTH_RELATIVE_FLOOR: 1
all-width gain label NONE: 5
all-width gain label SPACE_10_to_100: 2
all-width gain label SPACE_100: 1
all-width gain label SAS_100: 1
all-width gain label SAS_10_to_100: 1

## DIRECT_2D TOTAL SPACE_conservative budget
rows=12 distinct_case_ids=1 positive_width_intervals=12 zero_width_rows=0 finite_endpoints=1
positive_width_SAS_10_to_100=2 SAS_100=2 SPACE_10_to_100=3 SPACE_100=0
recorded finite bank budget span:   4.5209040E-06 to   1.9446000E-03
classification UNAVAILABLE: 2
classification NEITHER_100: 7
classification SAS_100: 2
classification BOTH_RELATIVE_FLOOR: 1
all-width gain label NONE: 5
all-width gain label SAS_10_to_100: 2
all-width gain label SPACE_10_to_100: 3
all-width gain label SAS_100: 2
