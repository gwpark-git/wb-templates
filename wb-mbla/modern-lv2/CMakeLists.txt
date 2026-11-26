cmake_minimum_required(VERSION 3.16)
project(wb_mbla_cpp VERSION 0.1.0 LANGUAGES CXX)
include(GNUInstallDirs)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
find_package(GSL 2.0 REQUIRED)
add_library(mbla_core STATIC
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
add_library(mbla_case_driver STATIC src/calculation_CP_profile.cpp)
target_link_libraries(mbla_case_driver PUBLIC mbla_core)
add_executable(mbla src/main.cpp)
target_link_libraries(mbla PRIVATE mbla_case_driver)

install(TARGETS mbla RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
# Collect the shared grant from a detached export, the parent package, or the collection.
set(mbla_new_license "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.NEW")
if(NOT EXISTS "${mbla_new_license}")
  set(mbla_new_license "${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE.NEW")
endif()
if(NOT EXISTS "${mbla_new_license}")
  set(mbla_new_license "${CMAKE_CURRENT_SOURCE_DIR}/../../LICENSE")
endif()
install(FILES "${mbla_new_license}" DESTINATION ${CMAKE_INSTALL_DOCDIR} RENAME LICENSE.NEW)
install(FILES LICENSE NOTICE THIRD_PARTY_NOTICES.org
  DESTINATION ${CMAKE_INSTALL_DOCDIR})

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
