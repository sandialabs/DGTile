set(APP_SOURCES
  hydro.cpp
  hydro_bounds.cpp
  hydro_dt.cpp
  hydro_face_integral.cpp
  hydro_fluxes.cpp
  hydro_ics.cpp
  hydro_volume_integral.cpp
  hydro_vtk.cpp
  lua.cpp
  main.cpp
  passive.cpp
  passive_ics.cpp
  passive_vtk.cpp
  run.cpp
  setup.cpp
)

set_source_files_properties(${APP_SOURCES} PROPERTIES LANGUAGE ${DGT_LANGUAGE})

add_executable(dgtapp ${APP_SOURCES})
target_compile_features(dgtapp PUBLIC cxx_std_17)
set_target_properties(dgtapp PROPERTIES ${DGT_LANGUAGE}_ARCHITECTURES "${DGT_ARCHITECTURES}")
target_link_libraries(dgtapp PRIVATE dgtile)
