# Fortran replay summary

Inputs: runs/paper_minimum_aggregate_02/cpu_main100_comparisons.csv and runs/paper_minimum_aggregate_02/cpu_main100_budgets.csv
Comparison CSV rows: 400
Distinct comparison case IDs: 100
Absolute-only evidence rows: 0
Budget CSV rows: 157
Distinct budget case IDs: 2

Counts describe exact saved rows and distinct case IDs within each group. Cases recur across families, modes, caps and budget directions; group counts must not be pooled as independent cases.
No population confidence intervals are inferred. Speed ranges describe recorded cost ratios. A cost winner can remain at the reporting floor; it is not a resolved precision claim.

## DIRECT_2D QUERY ACTUAL_ERROR cap=  1.0000000E-06
rows=100 distinct_case_ids=100
unmatched=0 unreached=0 floor_pairs=0 resolved_pairs=100 one_floor_pairs=0
SAS_cost_wins=100 SPACE_cost_wins=0 resolved_SAS_wins=100 resolved_SPACE_wins=0 ties_or_vetoed=0
separated_winners=100 rows_with_dominance_veto=0
SPACE_dominates_SAS_rows=0 SAS_dominates_SPACE_rows=99
SPACE/SAS cost ratio (selected pairs):   7.3357748E+00 to   1.1338897E+02
SAS winner speed ratio:   7.3357748E+00 to   1.1338897E+02
SPACE winner speed ratio: unavailable
status ACTUAL_ERROR_MATCH: 100

## DIRECT_2D QUERY TARGET cap=  1.0000000E-06
rows=100 distinct_case_ids=100
unmatched=0 unreached=0 floor_pairs=0 resolved_pairs=100 one_floor_pairs=0
SAS_cost_wins=100 SPACE_cost_wins=0 resolved_SAS_wins=100 resolved_SPACE_wins=0 ties_or_vetoed=0
separated_winners=100 rows_with_dominance_veto=0
SPACE_dominates_SAS_rows=0 SAS_dominates_SPACE_rows=0
SPACE/SAS cost ratio (selected pairs):   7.0426555E+00 to   1.1338897E+02
SAS winner speed ratio:   7.0426555E+00 to   1.1338897E+02
SPACE winner speed ratio: unavailable
status TARGET_BOTH_REACHED: 100

## DIRECT_2D TOTAL ACTUAL_ERROR cap=  1.0000000E-06
rows=100 distinct_case_ids=100
unmatched=0 unreached=0 floor_pairs=0 resolved_pairs=100 one_floor_pairs=0
SAS_cost_wins=100 SPACE_cost_wins=0 resolved_SAS_wins=100 resolved_SPACE_wins=0 ties_or_vetoed=0
separated_winners=100 rows_with_dominance_veto=0
SPACE_dominates_SAS_rows=0 SAS_dominates_SPACE_rows=99
SPACE/SAS cost ratio (selected pairs):   7.1863357E+00 to   1.1570470E+02
SAS winner speed ratio:   7.1863357E+00 to   1.1570470E+02
SPACE winner speed ratio: unavailable
status ACTUAL_ERROR_MATCH: 100

## DIRECT_2D TOTAL TARGET cap=  1.0000000E-06
rows=100 distinct_case_ids=100
unmatched=0 unreached=0 floor_pairs=0 resolved_pairs=100 one_floor_pairs=0
SAS_cost_wins=100 SPACE_cost_wins=0 resolved_SAS_wins=100 resolved_SPACE_wins=0 ties_or_vetoed=0
separated_winners=100 rows_with_dominance_veto=0
SPACE_dominates_SAS_rows=0 SAS_dominates_SPACE_rows=0
SPACE/SAS cost ratio (selected pairs):   7.1863357E+00 to   1.1570470E+02
SAS winner speed ratio:   7.1863357E+00 to   1.1570470E+02
SPACE winner speed ratio: unavailable
status TARGET_BOTH_REACHED: 100

Budget gains count only intervals with end > start: 10-to-100 means >=10 and <100; 100 means >=100. Finite terminal endpoints and zero-width rows are retained separately.

## DIRECT_2D QUERY nominal budget
rows=28 distinct_case_ids=2 positive_width_intervals=28 zero_width_rows=0 finite_endpoints=2
positive_width_SAS_10_to_100=0 SAS_100=14 SPACE_10_to_100=0 SPACE_100=0
recorded finite bank budget span:   9.4481132E-06 to   1.8843000E-02
classification UNAVAILABLE: 12
classification SAS_100: 14
classification BOTH_RELATIVE_FLOOR: 2
all-width gain label NONE: 14
all-width gain label SAS_100: 14

## DIRECT_2D QUERY SAS_conservative budget
rows=26 distinct_case_ids=2 positive_width_intervals=26 zero_width_rows=0 finite_endpoints=2
positive_width_SAS_10_to_100=0 SAS_100=12 SPACE_10_to_100=0 SPACE_100=0
recorded finite bank budget span:   9.4481132E-06 to   1.8843000E-02
classification UNAVAILABLE: 12
classification SAS_100: 12
classification BOTH_RELATIVE_FLOOR: 2
all-width gain label NONE: 14
all-width gain label SAS_100: 12

## DIRECT_2D QUERY SPACE_conservative budget
rows=24 distinct_case_ids=2 positive_width_intervals=24 zero_width_rows=0 finite_endpoints=2
positive_width_SAS_10_to_100=0 SAS_100=12 SPACE_10_to_100=0 SPACE_100=0
recorded finite bank budget span:   9.4481132E-06 to   1.8843000E-02
classification UNAVAILABLE: 10
classification SAS_100: 12
classification BOTH_RELATIVE_FLOOR: 2
all-width gain label NONE: 12
all-width gain label SAS_100: 12

## DIRECT_2D TOTAL nominal budget
rows=27 distinct_case_ids=2 positive_width_intervals=27 zero_width_rows=0 finite_endpoints=2
positive_width_SAS_10_to_100=0 SAS_100=13 SPACE_10_to_100=0 SPACE_100=0
recorded finite bank budget span:   9.8938272E-06 to   1.8203000E-02
classification UNAVAILABLE: 12
classification SAS_100: 13
classification BOTH_RELATIVE_FLOOR: 2
all-width gain label NONE: 14
all-width gain label SAS_100: 13

## DIRECT_2D TOTAL SAS_conservative budget
rows=28 distinct_case_ids=2 positive_width_intervals=28 zero_width_rows=0 finite_endpoints=2
positive_width_SAS_10_to_100=0 SAS_100=15 SPACE_10_to_100=0 SPACE_100=0
recorded finite bank budget span:   9.8938272E-06 to   1.8203000E-02
classification UNAVAILABLE: 11
classification SAS_100: 15
classification BOTH_RELATIVE_FLOOR: 2
all-width gain label NONE: 13
all-width gain label SAS_100: 15

## DIRECT_2D TOTAL SPACE_conservative budget
rows=24 distinct_case_ids=2 positive_width_intervals=24 zero_width_rows=0 finite_endpoints=2
positive_width_SAS_10_to_100=0 SAS_100=11 SPACE_10_to_100=0 SPACE_100=0
recorded finite bank budget span:   9.8938272E-06 to   1.8203000E-02
classification UNAVAILABLE: 11
classification SAS_100: 11
classification BOTH_RELATIVE_FLOOR: 2
all-width gain label NONE: 13
all-width gain label SAS_100: 11
