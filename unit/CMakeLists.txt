set(DG_UNIT_TEST_SOURCES
  grid.cpp
  spatial.cpp
  amr.cpp
  basis.cpp
  tree.cpp
  mesh.cpp
  file.cpp
  unit_tests.cpp
)

set_source_files_properties(${DG_UNIT_TEST_SOURCES} PROPERTIES LANGUAGE ${p3a_LANGUAGE})
add_executable(dgt-unit-tests ${DG_UNIT_TEST_SOURCES})
set_target_properties(dgt-unit-tests PROPERTIES ${p3a_LANGUAGE}_ARCHITECTURES "${p3a_ARCHITECTURES}")
target_link_libraries(dgt-unit-tests PRIVATE dgtile)
target_link_libraries(dgt-unit-tests PRIVATE GTest::gtest_main)
add_test(NAME unit-tests COMMAND dgt-unit-tests)
