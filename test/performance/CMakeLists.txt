set(PERFORMANCE_TEST_SOURCES
  packing.cpp
)

set_source_files_properties(${PERFORMANCE_TEST_SOURCES} PROPERTIES LANGUAGE ${DGT_LANGUAGE})

add_executable(dgt-performance-tests ${PERFORMANCE_TEST_SOURCES})
target_compile_features(dgt-performance-tests PUBLIC cxx_std_17)
set_target_properties(dgt-performance-tests PROPERTIES ${DGT_LANGUAGE}_ARCHITECTURES "${DGT_ARCHITECTURES}")
target_link_libraries(dgt-performance-tests PRIVATE dgtile)
