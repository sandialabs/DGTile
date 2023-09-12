set(INTEGRATION_TEST_SOURCES
  basis.cpp
  main.cpp
  when.cpp
)

set_source_files_properties(${INTEGRATION_TEST_SOURCES} PROPERTIES LANGUAGE ${DGT_LANGUAGE})

add_executable(dgt-integration-tests ${INTEGRATION_TEST_SOURCES})
target_compile_features(dgt-integration-tests PUBLIC cxx_std_17)
set_target_properties(dgt-integration-tests PROPERTIES ${DGT_LANGUAGE}_ARCHITECTURES "${DGT_ARCHITECTURES}")
target_link_libraries(dgt-integration-tests PRIVATE GTest::gtest)
target_link_libraries(dgt-integration-tests PRIVATE dgtile)

add_test(NAME integration_tests COMMAND "./dgt-integration-tests")
set_tests_properties(integration_tests PROPERTIES ENVIRONMENT "INTEGRATION_DATA_DIR=${CMAKE_CURRENT_SOURCE_DIR}")
