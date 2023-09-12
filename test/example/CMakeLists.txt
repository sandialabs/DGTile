set(EXAMPLE_SOURCES
  main.cpp
)

set_source_files_properties(${EXAMPLE_SOURCES} PROPERTIES LANGUAGE ${DGT_LANGUAGE})

add_executable(dgt-example ${EXAMPLE_SOURCES})
target_compile_features(dgt-example PUBLIC cxx_std_17)
set_target_properties(dgt-example PROPERTIES ${DGT_LANGUAGE}_ARCHITECTURES "${DGT_ARCHITECTURES}")
target_link_libraries(dgt-example PRIVATE dgtile)
