set(UNIT_TEST_SOURCES
  array.cpp
  base64.cpp
  block_info.cpp
  bounds.cpp
  box3.cpp
  cartesian.cpp
  defines.cpp
  dg.cpp
  dg_dim1.cpp
  dg_dim2.cpp
  dg_dim3.cpp
  field.cpp
  for_each.cpp
  ghosting.cpp
  grid3.cpp
  integrator.cpp
  interpolate.cpp
  mesh.cpp
  main.cpp
  partitioning.cpp
  print.cpp
  reduce.cpp
  stl.cpp
  subgrid3.cpp
  tree.cpp
  tree_amr.cpp
  tree_uniform.cpp
  vec.cpp
  vec3.cpp
  view.cpp
  when.cpp
)

if(DGT_ENABLE_CONDUIT)
  list(APPEND UNIT_TEST_SOURCES conduit.cpp)
endif()

set_source_files_properties(${UNIT_TEST_SOURCES} PROPERTIES LANGUAGE ${DGT_LANGUAGE})

add_executable(dgt-unit-tests ${UNIT_TEST_SOURCES})
target_compile_features(dgt-unit-tests PUBLIC cxx_std_17)
set_target_properties(dgt-unit-tests PROPERTIES ${DGT_LANGUAGE}_ARCHITECTURES "${DGT_ARCHITECTURES}")
target_link_libraries(dgt-unit-tests PRIVATE GTest::gtest)
target_link_libraries(dgt-unit-tests PRIVATE dgtile)

add_test(NAME unit_tests COMMAND "./dgt-unit-tests")
set_tests_properties(unit_tests PROPERTIES ENVIRONMENT "UNIT_DATA_DIR=${CMAKE_CURRENT_SOURCE_DIR}")
