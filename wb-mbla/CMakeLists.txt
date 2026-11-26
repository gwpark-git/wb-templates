cmake_minimum_required(VERSION 3.16)
project(wb_mbla_cpp VERSION 0.1.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
find_package(GSL 2.0 REQUIRED)
add_library(mbla_core
  src/dispersion_properties.cpp src/file_IO.cpp src/matrix.cpp
  src/osmotic_pressure_CS.cpp src/parameters.cpp src/ported_lib.cpp
  src/transport_properties_SPHS.cpp src/transport_properties_SID.cpp
  src/transport_properties_INTERPOLATE.cpp src/boost_matrix.cpp)
target_include_directories(mbla_core PUBLIC src)
target_link_libraries(mbla_core PUBLIC GSL::gsl GSL::gslcblas)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  # Bounds-only edits changed default cake results on AArch64. A common
  # no-contraction policy recovered source parity; see the retained diagnosis.
  target_compile_options(mbla_core PUBLIC -ffp-contract=off)
else()
  message(WARNING "Unverified FP contraction policy: establish target-specific cake references")
endif()
add_executable(mbla src/calculation_CP_profile.cpp)
target_link_libraries(mbla PRIVATE mbla_core)
option(MBLA_SANITIZERS "Enable GCC/Clang address and undefined behavior sanitizers" OFF)
if(MBLA_SANITIZERS)
  target_compile_options(mbla_core PUBLIC -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_options(mbla_core PUBLIC -fsanitize=address,undefined)
endif()
include(CTest)
if(BUILD_TESTING)
  add_executable(mbla_units tests/unit.cpp)
  target_link_libraries(mbla_units PRIVATE mbla_core)
  add_test(NAME mbla_units COMMAND mbla_units)
  find_package(Python3 COMPONENTS Interpreter REQUIRED)
  add_test(NAME mbla_regression COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression.py --exe $<TARGET_FILE:mbla>
    --work ${CMAKE_CURRENT_BINARY_DIR}/test-runs)
  set_tests_properties(mbla_regression PROPERTIES TIMEOUT 90)
  foreach(contract IN ITEMS input output)
    add_test(NAME mbla_${contract}_contract COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/boundaries.py --exe $<TARGET_FILE:mbla>
      --work ${CMAKE_CURRENT_BINARY_DIR}/boundary-runs/${contract}
      --section ${contract})
    set_tests_properties(mbla_${contract}_contract PROPERTIES TIMEOUT 45)
  endforeach()
endif()
