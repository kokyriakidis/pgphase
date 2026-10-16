.DEFAULT_GOAL := all

CXX ?= g++
CC ?= gcc
THIRD_PARTY ?= $(abspath third_party)
WFA2_ROOT ?= $(THIRD_PARTY)/WFA2-lib
ABPOA_ROOT ?= $(THIRD_PARTY)/abPOA
EDLIB_ROOT ?= $(THIRD_PARTY)/edlib/edlib
# Rust is needed for pgcat (make pgcat). Prefer the rustup-managed toolchain
# over any system-installed cargo.
CARGO ?= $(firstword $(wildcard $(HOME)/.cargo/bin/cargo) cargo)

CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -MMD -MP
WFA_CPPFLAGS = -I$(WFA2_ROOT)
AB_CPPFLAGS = -I$(ABPOA_ROOT)/include
EDLIB_CPPFLAGS = -I$(EDLIB_ROOT)/include
C_CFLAGS = -O3 -Wall -Wextra
ALIGN_CPPFLAGS = $(WFA_CPPFLAGS) $(AB_CPPFLAGS) $(EDLIB_CPPFLAGS)

WFA2_LIB = $(WFA2_ROOT)/lib/libwfa.a
ABPOA_LIB = $(ABPOA_ROOT)/lib/libabpoa.a

SOURCES_CXX = src/main.cpp \
	src/allele_identity.cpp \
	src/collect_pipeline.cpp \
	src/pgcat_command.cpp \
	src/graph_collect.cpp \
	src/bam_digar.cpp \
	src/noise_filter.cpp \
	src/collect_var.cpp \
	src/collect_phase.cpp \
	src/collect_phase_noisy.cpp \
	src/align.cpp \
	src/collect_bam_output.cpp \
	src/collect_output.cpp \
	src/graph_bam_adapter.cpp \
	src/graph_sites.cpp \
	src/graph_query.cpp \
	src/union_phase.cpp \
	src/union_inject.cpp \
	src/union_em.cpp \
	src/union_windows.cpp \
	src/union_reads.cpp
SOURCES_C = src/sdust.c src/cgranges.c src/kalloc.c

OBJS = $(SOURCES_CXX:.cpp=.o) $(SOURCES_C:.c=.o)
UNION_OBJS = src/union_phase.o src/union_inject.o src/union_em.o src/union_windows.o src/union_reads.o
EDLIB_OBJ = src/edlib.o
OBJS += $(EDLIB_OBJ)

LDFLAGS ?= -lhts -lzstd -lm -lz -lpthread

-include $(patsubst %.cpp,%.d,$(SOURCES_CXX))

.PHONY: all clean check unit-tests predicate-tests parity-tests pgcat-tests upstream-parity-tests gap-dev-check third-party-libs portable-bundle release release-strict

all: pgphase

check: pgphase
	bash scripts/validate_collect_gates.sh

predicate-tests: test_phase_predicates
	./test_phase_predicates

# In-memory fixtures are the edit loop; owning replays are an explicit check.
.PHONY: gap-dev-check
gap-dev-check: test_phase_predicates test_graph_bam_adapter test_allele_identity
	./test_allele_identity
	./test_phase_predicates $(if $(PREDICATE),"$(PREDICATE)",)
	./test_graph_bam_adapter

parity-tests: test_port_parity
	./test_port_parity

# Requires a local longcallD checkout; compiles its original C source so the
# expected results cannot drift with a hand-copied reference implementation.
LONGCALLD_ROOT ?= $(abspath ../longcallD)
UPSTREAM_PHASE_DIR = build/upstream_phase
UPSTREAM_PHASE_INCLUDES = -I$(LONGCALLD_ROOT)/src -I$(LONGCALLD_ROOT)/htslib -I$(LONGCALLD_ROOT)/edlib/include -I$(LONGCALLD_ROOT)/abPOA/include -I$(LONGCALLD_ROOT)/WFA2-lib
-include $(UPSTREAM_PHASE_DIR)/assign_hap.d $(UPSTREAM_PHASE_DIR)/collect_var.d $(UPSTREAM_PHASE_DIR)/bridge.d

upstream-parity-tests: test_upstream_phase
	./test_upstream_phase

$(UPSTREAM_PHASE_DIR):
	mkdir -p $@

$(UPSTREAM_PHASE_DIR)/assign_hap.o: $(LONGCALLD_ROOT)/src/assign_hap.c | $(UPSTREAM_PHASE_DIR)
	$(CC) $(C_CFLAGS) -MMD -MP -Wno-unused-parameter -Wno-unused-function -Wno-misleading-indentation -ffunction-sections -fdata-sections $(UPSTREAM_PHASE_INCLUDES) -c $< -o $@

$(UPSTREAM_PHASE_DIR)/collect_var.o: $(LONGCALLD_ROOT)/src/collect_var.c | $(UPSTREAM_PHASE_DIR)
	$(CC) $(C_CFLAGS) -MMD -MP -Wno-unused-parameter -Wno-unused-function -Wno-misleading-indentation -Wno-absolute-value -ffunction-sections -fdata-sections $(UPSTREAM_PHASE_INCLUDES) -c $< -o $@

$(UPSTREAM_PHASE_DIR)/bridge.o: src/test_upstream_phase_bridge.c src/test_upstream_phase_bridge.h | $(UPSTREAM_PHASE_DIR)
	$(CC) $(C_CFLAGS) -MMD -MP -Wno-unused-function $(UPSTREAM_PHASE_INCLUDES) -Isrc -c $< -o $@

unit-tests: test_graph_sites test_graph_bam_adapter test_noise_filter test_allele_identity test_union_phase
	./test_allele_identity
	./test_graph_sites
	./test_graph_bam_adapter
	./test_noise_filter
	./test_union_phase

portable-bundle: pgphase
	bash scripts/make_portable_bundle.sh

release: pgphase
	bash scripts/make_release_bundle.sh

release-strict: pgphase
	RUN_CHECKS=1 bash scripts/make_release_bundle.sh

src/align.o: src/align.cpp
	$(CXX) $(CXXFLAGS) $(ALIGN_CPPFLAGS) -c $< -o $@

src/graph_collect.o: src/graph_collect.cpp
	$(CXX) $(CXXFLAGS) $(EDLIB_CPPFLAGS) -c $< -o $@

$(UNION_OBJS): src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) $(EDLIB_CPPFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(C_CFLAGS) -c $< -o $@

third-party-libs: $(WFA2_LIB) $(ABPOA_LIB)

$(WFA2_LIB):
	$(MAKE) -C "$(WFA2_ROOT)" setup
	$(MAKE) -C "$(WFA2_ROOT)" lib_wfa

$(ABPOA_LIB):
	$(MAKE) -C "$(ABPOA_ROOT)" libabpoa

$(EDLIB_OBJ): $(EDLIB_ROOT)/src/edlib.cpp
	$(CXX) $(CXXFLAGS) $(EDLIB_CPPFLAGS) -c $< -o $@

pgphase: $(OBJS) $(WFA2_LIB) $(ABPOA_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(WFA2_LIB) $(ABPOA_LIB) $(LDFLAGS)

# pgcat builds the site catalogs (pgphase build-snarl-catalog runs it); it is
# placed beside pgphase, where build-snarl-catalog looks first.
PGCAT_SOURCES = tools/pgcat/Cargo.toml tools/pgcat/Cargo.lock $(wildcard tools/pgcat/src/*.rs)
pgcat: $(PGCAT_SOURCES)
	cd tools/pgcat && $(CARGO) build --release --locked
	cp tools/pgcat/target/release/pgcat $@

pgcat-tests:
	cd tools/pgcat && $(CARGO) test --release --locked


test_graph_sites: src/test_graph_sites.cpp src/graph_sites.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

test_allele_identity: src/test_allele_identity.cpp src/allele_identity.o
	$(CXX) $(CXXFLAGS) -o $@ $^

test_graph_bam_adapter: src/test_graph_bam_adapter.cpp src/graph_bam_adapter.o src/allele_identity.o src/noise_filter.o src/graph_sites.o src/graph_query.o src/collect_phase.o src/collect_output.o src/cgranges.o src/kalloc.o src/sdust.o src/collect_var.o src/collect_phase_noisy.o src/align.o $(EDLIB_OBJ) $(WFA2_LIB) $(ABPOA_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

test_union_phase: src/test_union_phase.cpp $(UNION_OBJS) src/graph_bam_adapter.o src/allele_identity.o src/noise_filter.o src/graph_sites.o src/graph_query.o src/collect_phase.o src/collect_output.o src/cgranges.o src/kalloc.o src/sdust.o src/collect_var.o src/collect_phase_noisy.o src/align.o $(EDLIB_OBJ) $(WFA2_LIB) $(ABPOA_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

test_noise_filter: src/test_noise_filter.cpp src/graph_bam_adapter.o src/allele_identity.o src/noise_filter.o src/graph_sites.o src/graph_query.o src/collect_phase.o src/collect_output.o src/cgranges.o src/kalloc.o src/sdust.o src/collect_var.o src/collect_phase_noisy.o src/align.o $(EDLIB_OBJ) $(WFA2_LIB) $(ABPOA_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Catch2 v2 is a single vendored header; -O1 keeps its compile time tolerable.
test_phase_predicates: src/test_phase_predicates.cpp third_party/catch2/catch.hpp src/collect_phase.o src/collect_phase_noisy.o src/collect_output.o src/collect_var.o src/noise_filter.o src/align.o src/cgranges.o src/kalloc.o src/sdust.o $(EDLIB_OBJ) $(WFA2_LIB) $(ABPOA_LIB)
	$(CXX) -O1 -std=c++17 -Wall -Wextra -Isrc -I. -o $@ $< src/collect_phase.o src/collect_phase_noisy.o src/collect_output.o src/collect_var.o src/noise_filter.o src/align.o src/cgranges.o src/kalloc.o src/sdust.o $(EDLIB_OBJ) $(WFA2_LIB) $(ABPOA_LIB) $(LDFLAGS)


clean:
	rm -f pgphase pgcat test_union_phase test_phase_predicates test_port_parity test_upstream_phase test_graph_sites test_graph_bam_adapter test_noise_filter test_allele_identity src/*.o src/*.d
	rm -rf $(UPSTREAM_PHASE_DIR) tools/pgcat/target

test_upstream_phase: src/test_upstream_phase.cpp src/test_upstream_phase_bridge.h $(UPSTREAM_PHASE_DIR)/bridge.o $(UPSTREAM_PHASE_DIR)/assign_hap.o $(UPSTREAM_PHASE_DIR)/collect_var.o src/collect_phase.o src/collect_phase_noisy.o src/collect_output.o src/collect_var.o src/noise_filter.o src/align.o src/cgranges.o src/kalloc.o src/sdust.o $(EDLIB_OBJ) $(WFA2_LIB) $(ABPOA_LIB)
	$(CXX) -O1 -std=c++17 -Wall -Wextra -Isrc -I. -o $@ $< $(UPSTREAM_PHASE_DIR)/bridge.o $(UPSTREAM_PHASE_DIR)/assign_hap.o $(UPSTREAM_PHASE_DIR)/collect_var.o src/collect_phase.o src/collect_phase_noisy.o src/collect_output.o src/collect_var.o src/noise_filter.o src/align.o src/cgranges.o src/kalloc.o src/sdust.o $(EDLIB_OBJ) $(WFA2_LIB) $(ABPOA_LIB) -Wl,--gc-sections $(LDFLAGS)

test_port_parity: src/test_port_parity.cpp third_party/catch2/catch.hpp src/collect_phase.o src/collect_phase_noisy.o src/collect_output.o src/collect_var.o src/noise_filter.o src/align.o src/cgranges.o src/kalloc.o src/sdust.o $(EDLIB_OBJ) $(WFA2_LIB) $(ABPOA_LIB)
	$(CXX) -O1 -std=c++17 -Wall -Wextra -Isrc -I. -o $@ $< src/collect_phase.o src/collect_phase_noisy.o src/collect_output.o src/collect_var.o src/noise_filter.o src/align.o src/cgranges.o src/kalloc.o src/sdust.o $(EDLIB_OBJ) $(WFA2_LIB) $(ABPOA_LIB) $(LDFLAGS)
