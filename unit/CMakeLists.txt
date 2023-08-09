set(UNIT_TEST_SOURCES
  base64.cpp
  bounds.cpp
  box3.cpp
  cartesian.cpp
  defines.cpp
  dg.cpp
  for_each.cpp
  grid3.cpp
  main.cpp
  subgrid3.cpp
  vec.cpp
  vec3.cpp
  when.cpp
)

set_source_files_properties(${UNIT_TEST_SOURCES} PROPERTIES LANGUAGE ${DGT_LANGUAGE})

add_executable(dgt-unit-tests ${UNIT_TEST_SOURCES})
target_compile_features(dgt-unit-tests PUBLIC cxx_std_17)
set_target_properties(dgt-unit-tests PROPERTIES ${DGT_LANGUAGE}_ARCHITECTURES "${DGT_ARCHITECTURES}")
target_link_libraries(dgt-unit-tests PRIVATE GTest::gtest)
target_link_libraries(dgt-unit-tests PRIVATE dgtile)

add_test(NAME unit_tests COMMAND "./dgt-unit-tests")
set_tests_properties(unit_tests PROPERTIES ENVIRONMENT "DATA_DIR=${CMAKE_CURRENT_SOURCE_DIR}")
