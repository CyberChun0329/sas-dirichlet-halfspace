# Compile only. No target executes a scientific, test, timing, or reducer binary.
FC = gfortran
BASE = -std=f2008 -Wall -Wextra -Wno-compare-reals -Wimplicit-interface -fimplicit-none -ffree-line-length-none
RELEASE = $(BASE) -O3 -march=native
DEBUG = $(BASE) -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
MPFR_PREFIX ?= /opt/homebrew/opt/mpfr
GMP_PREFIX ?= /opt/homebrew/opt/gmp
MPFR_LIBS ?= -L$(MPFR_PREFIX)/lib -lmpfr -L$(GMP_PREFIX)/lib -lgmp
ACCELERATE ?= -framework Accelerate

CORE = vendor/core/temperature_common.f90 vendor/core/prepared_sas.f90 vendor/core/prepared_space.f90
MATH = $(CORE) vendor/performance/replay_types.f90 vendor/performance/gaussian_axes.f90 vendor/performance/sas_general.f90 vendor/performance/spatial_routes.f90 vendor/performance/source_projection.f90 vendor/performance/replay_measure.f90 vendor/performance/replay_select.f90 vendor/performance/svd_factor_routes.f90 vendor/performance/replay_case_bank.f90 vendor/performance/useful_screening.f90
MP = src/validation/mp_arithmetic.f90 src/validation/mp_moments.f90
VALIDATION = vendor/core/three_dimensional_analysis.f90 vendor/core/temperature_common.f90 vendor/performance/replay_types.f90 vendor/performance/gaussian_axes.f90 src/validation/validation_sas.f90 $(MP) src/validation/spatial_reference.f90
ENGINE = $(MATH) src/unified_policy.f90 src/performance/unified_engine.f90 src/performance/paired_measure.f90
POLY = $(ENGINE) src/performance/unified_poly_reference.f90
ORIGINAL_CONTROLS = $(CORE) vendor/performance/replay_types.f90 vendor/performance/gaussian_axes.f90 vendor/performance/sas_general.f90 vendor/performance/spatial_routes.f90 vendor/performance/source_projection.f90

PAPER_BINS = validation_driver precision_driver periodic_driver summarize_validation stage_c_reduce gaussian_coeff_check gaussian_original_controls targeted_driver replay_summary
CORE_TESTS = temperature_unit test_sas_general test_axis_recurrence test_spatial_routes test_source_projection test_replay_select test_replay_measure test_svd_factor_routes test_case_bank
TEST_BINS = $(CORE_TESTS) validation_smoke test_policy test_unified_engine test_paired_measure test_unified_reference

.DEFAULT_GOAL := all
.PHONY: all release debug build-tests build-reference-generators
all: release
release: $(addprefix build/release/,$(PAPER_BINS))
debug: $(addprefix build/debug/,$(PAPER_BINS))
build-tests: $(addprefix build/debug/,$(TEST_BINS)) $(addprefix build/release/,$(TEST_BINS))
# Reference generation is an opt-in compilation target.
build-reference-generators: build/release/unified_poly_driver

# Separate module directories prevent races with make -j2. All build products
# (including debug symbols) stay below build; copied science sources stay intact.
define executable
build/release/$(1): $(2) Makefile
	mkdir -p build/release build/mod/release/$(1)
	$(FC) $(RELEASE) -Jbuild/mod/release/$(1) -Ibuild/mod/release/$(1) $(2) $(3) -o $$@
build/debug/$(1): $(2) Makefile
	mkdir -p build/debug build/mod/debug/$(1)
	$(FC) $(DEBUG) -Jbuild/mod/debug/$(1) -Ibuild/mod/debug/$(1) $(2) $(3) -o $$@
endef

$(eval $(call executable,validation_driver,$(VALIDATION) src/validation/validation_driver.f90,$(MPFR_LIBS)))
$(eval $(call executable,precision_driver,$(MP) src/validation/precision_driver.f90,$(MPFR_LIBS)))
$(eval $(call executable,periodic_driver,$(VALIDATION) src/validation/periodic_validation.f90 src/validation/periodic_driver.f90,$(MPFR_LIBS)))
$(eval $(call executable,summarize_validation,src/validation/summarize_validation.f90))
$(eval $(call executable,stage_c_reduce,stage_c_support/stage_c_reduce.f90))
$(eval $(call executable,gaussian_coeff_check,$(MP) stage_c_support/gaussian_coeff_check.f90,$(MPFR_LIBS)))
$(eval $(call executable,gaussian_original_controls,$(ORIGINAL_CONTROLS) stage_c_support/gaussian_original_controls.f90,$(ACCELERATE)))
$(eval $(call executable,targeted_driver,$(POLY) paper_minimum/targeted_driver.f90,$(ACCELERATE)))
$(eval $(call executable,replay_summary,paper_minimum/replay_summary.f90))
$(eval $(call executable,unified_poly_driver,$(POLY) src/performance/unified_poly_driver.f90,$(ACCELERATE)))

$(foreach name,$(CORE_TESTS),$(eval $(call executable,$(name),$(MATH) tests/core/$(name).f90,$(ACCELERATE))))
$(eval $(call executable,validation_smoke,$(VALIDATION) tests/validation/validation_smoke.f90,$(MPFR_LIBS)))
$(eval $(call executable,test_policy,src/unified_policy.f90 tests/test_policy.f90))
$(foreach name,test_unified_engine test_paired_measure,$(eval $(call executable,$(name),$(ENGINE) tests/performance/$(name).f90,$(ACCELERATE))))
$(eval $(call executable,test_unified_reference,$(POLY) tests/performance/test_unified_reference.f90,$(ACCELERATE)))
