# Fortran replay summary

Inputs: runs/paper_minimum_aggregate_02/wall_budget3_comparisons.csv and runs/paper_minimum_aggregate_02/wall_budget3_budgets.csv
Comparison CSV rows: 12
Distinct comparison case IDs: 3
Absolute-only evidence rows: 0
Budget CSV rows: 223
Distinct budget case IDs: 3

Counts describe exact saved rows and distinct case IDs within each group. Cases recur across families, modes, caps and budget directions; group counts must not be pooled as independent cases.
No population confidence intervals are inferred. Speed ranges describe recorded cost ratios. A cost winner can remain at the reporting floor; it is not a resolved precision claim.

## DIRECT_2D QUERY ACTUAL_ERROR cap=  1.0000000E-06
rows=3 distinct_case_ids=3
unmatched=0 unreached=0 floor_pairs=0 resolved_pairs=3 one_floor_pairs=0
SAS_cost_wins=2 SPACE_cost_wins=1 resolved_SAS_wins=2 resolved_SPACE_wins=1 ties_or_vetoed=0
separated_winners=3 rows_with_dominance_veto=0
SPACE_dominates_SAS_rows=0 SAS_dominates_SPACE_rows=2
SPACE/SAS cost ratio (selected pairs):   7.7300780E-01 to   8.5779970E+01
SAS winner speed ratio:   1.1541494E+01 to   8.5779970E+01
SPACE winner speed ratio:   1.2936480E+00 to   1.2936480E+00
status ACTUAL_ERROR_MATCH: 3

## DIRECT_2D QUERY TARGET cap=  1.0000000E-06
rows=3 distinct_case_ids=3
unmatched=0 unreached=0 floor_pairs=0 resolved_pairs=3 one_floor_pairs=0
SAS_cost_wins=2 SPACE_cost_wins=1 resolved_SAS_wins=2 resolved_SPACE_wins=1 ties_or_vetoed=0
separated_winners=3 rows_with_dominance_veto=0
SPACE_dominates_SAS_rows=0 SAS_dominates_SPACE_rows=0
SPACE/SAS cost ratio (selected pairs):   7.7300780E-01 to   8.5779970E+01
SAS winner speed ratio:   1.1541494E+01 to   8.5779970E+01
SPACE winner speed ratio:   1.2936480E+00 to   1.2936480E+00
status TARGET_BOTH_REACHED: 3

## DIRECT_2D TOTAL ACTUAL_ERROR cap=  1.0000000E-06
rows=3 distinct_case_ids=3
unmatched=0 unreached=0 floor_pairs=0 resolved_pairs=3 one_floor_pairs=0
SAS_cost_wins=2 SPACE_cost_wins=1 resolved_SAS_wins=2 resolved_SPACE_wins=1 ties_or_vetoed=0
separated_winners=3 rows_with_dominance_veto=0
SPACE_dominates_SAS_rows=0 SAS_dominates_SPACE_rows=2
SPACE/SAS cost ratio (selected pairs):   8.2515933E-01 to   9.5156608E+01
SAS winner speed ratio:   1.0930176E+01 to   9.5156608E+01
SPACE winner speed ratio:   1.2118872E+00 to   1.2118872E+00
status ACTUAL_ERROR_MATCH: 3

## DIRECT_2D TOTAL TARGET cap=  1.0000000E-06
rows=3 distinct_case_ids=3
unmatched=0 unreached=0 floor_pairs=0 resolved_pairs=3 one_floor_pairs=0
SAS_cost_wins=2 SPACE_cost_wins=1 resolved_SAS_wins=2 resolved_SPACE_wins=1 ties_or_vetoed=0
separated_winners=3 rows_with_dominance_veto=0
SPACE_dominates_SAS_rows=0 SAS_dominates_SPACE_rows=0
SPACE/SAS cost ratio (selected pairs):   8.2515933E-01 to   9.5156608E+01
SAS winner speed ratio:   1.0930176E+01 to   9.5156608E+01
SPACE winner speed ratio:   1.2118872E+00 to   1.2118872E+00
status TARGET_BOTH_REACHED: 3

Budget gains count only intervals with end > start: 10-to-100 means >=10 and <100; 100 means >=100. Finite terminal endpoints and zero-width rows are retained separately.

## DIRECT_2D QUERY nominal budget
rows=39 distinct_case_ids=3 positive_width_intervals=39 zero_width_rows=0 finite_endpoints=3
positive_width_SAS_10_to_100=1 SAS_100=14 SPACE_10_to_100=3 SPACE_100=1
recorded finite bank budget span:   3.9564787E-06 to   2.2138000E-02
classification UNAVAILABLE: 15
classification SAS_100: 14
classification BOTH_RELATIVE_FLOOR: 3
classification NEITHER_100: 6
classification SPACE_100: 1
all-width gain label NONE: 20
all-width gain label SAS_100: 14
all-width gain label SPACE_10_to_100: 3
all-width gain label SPACE_100: 1
all-width gain label SAS_10_to_100: 1

## DIRECT_2D QUERY SAS_conservative budget
rows=37 distinct_case_ids=3 positive_width_intervals=37 zero_width_rows=0 finite_endpoints=3
positive_width_SAS_10_to_100=1 SAS_100=14 SPACE_10_to_100=3 SPACE_100=2
recorded finite bank budget span:   3.9564787E-06 to   2.2138000E-02
classification UNAVAILABLE: 13
classification SAS_100: 14
classification BOTH_RELATIVE_FLOOR: 3
classification NEITHER_100: 5
classification SPACE_100: 2
all-width gain label NONE: 17
all-width gain label SAS_100: 14
all-width gain label SPACE_10_to_100: 3
all-width gain label SPACE_100: 2
all-width gain label SAS_10_to_100: 1

## DIRECT_2D QUERY SPACE_conservative budget
rows=37 distinct_case_ids=3 positive_width_intervals=37 zero_width_rows=0 finite_endpoints=3
positive_width_SAS_10_to_100=1 SAS_100=14 SPACE_10_to_100=3 SPACE_100=0
recorded finite bank budget span:   3.9564787E-06 to   2.2138000E-02
classification UNAVAILABLE: 14
classification SAS_100: 14
classification BOTH_RELATIVE_FLOOR: 3
classification NEITHER_100: 6
all-width gain label NONE: 19
all-width gain label SAS_100: 14
all-width gain label SPACE_10_to_100: 3
all-width gain label SAS_10_to_100: 1

## DIRECT_2D TOTAL nominal budget
rows=39 distinct_case_ids=3 positive_width_intervals=39 zero_width_rows=0 finite_endpoints=3
positive_width_SAS_10_to_100=1 SAS_100=14 SPACE_10_to_100=3 SPACE_100=1
recorded finite bank budget span:   4.5209040E-06 to   1.8821000E-02
classification UNAVAILABLE: 15
classification SAS_100: 14
classification BOTH_RELATIVE_FLOOR: 3
classification NEITHER_100: 6
classification SPACE_100: 1
all-width gain label NONE: 20
all-width gain label SAS_100: 14
all-width gain label SPACE_10_to_100: 3
all-width gain label SPACE_100: 1
all-width gain label SAS_10_to_100: 1

## DIRECT_2D TOTAL SAS_conservative budget
rows=35 distinct_case_ids=3 positive_width_intervals=35 zero_width_rows=0 finite_endpoints=3
positive_width_SAS_10_to_100=1 SAS_100=15 SPACE_10_to_100=2 SPACE_100=1
recorded finite bank budget span:   4.5209040E-06 to   1.8821000E-02
classification UNAVAILABLE: 11
classification SAS_100: 15
classification BOTH_RELATIVE_FLOOR: 3
classification NEITHER_100: 5
classification SPACE_100: 1
all-width gain label NONE: 16
all-width gain label SAS_100: 15
all-width gain label SPACE_10_to_100: 2
all-width gain label SPACE_100: 1
all-width gain label SAS_10_to_100: 1

## DIRECT_2D TOTAL SPACE_conservative budget
rows=36 distinct_case_ids=3 positive_width_intervals=36 zero_width_rows=0 finite_endpoints=3
positive_width_SAS_10_to_100=2 SAS_100=13 SPACE_10_to_100=3 SPACE_100=0
recorded finite bank budget span:   4.5209040E-06 to   1.8821000E-02
classification UNAVAILABLE: 13
classification SAS_100: 13
classification BOTH_RELATIVE_FLOOR: 3
classification NEITHER_100: 7
all-width gain label NONE: 18
all-width gain label SAS_100: 13
all-width gain label SAS_10_to_100: 2
all-width gain label SPACE_10_to_100: 3
