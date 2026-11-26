cmake_minimum_required(VERSION 3.16)
project(wb_mbla_level3 LANGUAGES CXX)
include(GNUInstallDirs)
include(CTest)

# Reuse the existing case driver and numerical solver from the full mBLA template.
add_subdirectory(../modern-lv2 level2)
add_executable(mbla-lv3-unity examples/unity_provider.cpp)
target_compile_features(mbla-lv3-unity PRIVATE cxx_std_17)
target_link_libraries(mbla-lv3-unity PRIVATE mbla_case_driver)
install(TARGETS mbla-lv3-unity RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

if(BUILD_TESTING)
  find_package(Python3 COMPONENTS Interpreter REQUIRED)
  add_test(NAME mbla_lv3_unity_contract COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/compare_unity.py
    $<TARGET_FILE:mbla> $<TARGET_FILE:mbla-lv3-unity>
    ${CMAKE_CURRENT_SOURCE_DIR}/examples/transport-unity.inp)
  set_tests_properties(mbla_lv3_unity_contract PROPERTIES TIMEOUT 45)
endif()
