CXX ?= g++
PYTHON ?= python3
PKG_CONFIG ?= pkg-config
GSL_CFLAGS ?= $(shell $(PKG_CONFIG) --cflags gsl)
GSL_LIBS ?= $(shell $(PKG_CONFIG) --libs gsl)
OPENMP_FLAGS ?= -fopenmp

# Select a different directory when changing compilers or build flags.
BUILD_DIR ?= build/$(notdir $(CXX))
RUN_DIR ?= runs/regression
EXECUTABLE := $(BUILD_DIR)/stochastic_simulation
GIT_VERSION := $(shell if test -e .git; then git describe --dirty --always --tags 2>/dev/null || echo unknown; else echo 0.1.0-ai-modulated; fi)
CPPFLAGS += $(GSL_CFLAGS) -DVERSION=\"$(GIT_VERSION)\"
CXXFLAGS ?= -O2 -g
COMPILE_FLAGS := -std=c++17 -Wall -Wextra $(OPENMP_FLAGS)
LDLIBS += $(GSL_LIBS)

SOURCES := src/stochastic_simulation.cpp src/brownian.cpp src/dumbbell_model.cpp \
           src/repulsive_brownian.cpp src/stochastic_HEUR_flowers.cpp \
           lib/association.cpp lib/connectivity.cpp lib/geometry.cpp \
           lib/handle_association.cpp lib/matrix.cpp lib/parallel.cpp \
           lib/potential.cpp lib/random.cpp lib/file_IO.cpp \
           lib/time_evolution.cpp lib/trajectory.cpp lib/cell_list.cpp lib/analysis.cpp
OBJECTS := $(addprefix $(BUILD_DIR)/,$(SOURCES:.cpp=.o))
NUMERIC_TEST := $(BUILD_DIR)/test_numerics
TOPOLOGY_TEST := $(BUILD_DIR)/test_legacy_topology
PHYSICS_TEST := $(BUILD_DIR)/test_legacy_physics
ANALYSIS_TEST := $(BUILD_DIR)/test_legacy_analysis
ANALYSIS_TOOLS := $(BUILD_DIR)/correlation_function_individual_chain $(BUILD_DIR)/measure_stress_autocorrelation

.PHONY: all analysis-tools opt debug scope debug_scope_opt test test-unit check check-unit check-regression check-sanitize clean
all: $(EXECUTABLE) analysis-tools
analysis-tools: $(ANALYSIS_TOOLS)
opt scope: all
debug debug_scope_opt:
	$(MAKE) BUILD_DIR="$(BUILD_DIR)/debug" CXXFLAGS="-O0 -g3" all

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OPENMP_FLAGS) $^ $(LDLIBS) -o $@

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(COMPILE_FLAGS) -MMD -MP -c $< -o $@

$(NUMERIC_TEST): tests/test_numerics.cpp $(BUILD_DIR)/lib/matrix.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(COMPILE_FLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(TOPOLOGY_TEST): tests/test_legacy_topology.cpp $(filter-out $(BUILD_DIR)/src/stochastic_simulation.o,$(OBJECTS))
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(COMPILE_FLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(PHYSICS_TEST): tests/test_legacy_physics.cpp $(filter-out $(BUILD_DIR)/src/stochastic_simulation.o,$(OBJECTS))
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(COMPILE_FLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(ANALYSIS_TOOLS): $(BUILD_DIR)/%: src/%.cpp $(filter-out $(BUILD_DIR)/src/stochastic_simulation.o,$(OBJECTS))
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(COMPILE_FLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(ANALYSIS_TEST): tests/test_legacy_analysis.cpp $(BUILD_DIR)/lib/analysis.o $(BUILD_DIR)/lib/matrix.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(COMPILE_FLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

test-unit: $(NUMERIC_TEST) $(TOPOLOGY_TEST) $(PHYSICS_TEST) $(ANALYSIS_TEST)
	$(NUMERIC_TEST)
	$(TOPOLOGY_TEST)
	$(PHYSICS_TEST)
	$(ANALYSIS_TEST)

test: $(EXECUTABLE) test-unit
	$(PYTHON) tests/regression.py --executable $(EXECUTABLE) --run-dir $(RUN_DIR) $(if $(COMPARE_EXECUTABLE),--compare-executable $(COMPARE_EXECUTABLE))

# Unified, isolated validation: no package installation or baseline regeneration.
check:
	$(PYTHON) tools/validate.py --compiler $(CXX) --profile full
check-unit:
	$(PYTHON) tools/validate.py --compiler $(CXX) --profile unit
check-regression:
	$(PYTHON) tools/validate.py --compiler $(CXX) --profile regression
check-sanitize:
	$(PYTHON) tools/validate.py --compiler $(CXX) --profile sanitizer

clean:
	$(RM) $(OBJECTS) $(OBJECTS:.o=.d) $(EXECUTABLE) $(NUMERIC_TEST) $(TOPOLOGY_TEST) $(PHYSICS_TEST) $(ANALYSIS_TEST) $(ANALYSIS_TOOLS)

-include $(OBJECTS:.o=.d)
