set(APP_SOURCES
  lua.cpp
  main.cpp
  run.cpp
  setup.cpp
)

set_source_files_properties(${APP_SOURCES} PROPERTIES LANGUAGE ${DGT_LANGUAGE})

add_executable(dgtapp ${APP_SOURCES})
target_compile_features(dgtapp PUBLIC cxx_std_17)
set_target_properties(dgtapp PROPERTIES ${DGT_LANGUAGE}_ARCHITECTURES "${DGT_ARCHITECTURES}")
target_link_libraries(dgtapp PRIVATE dgtile)
