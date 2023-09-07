set(DGT_HEADERS
  dgt_base64.hpp
  dgt_bounds.hpp
  dgt_box3.hpp
  dgt_cartesian.hpp
  dgt_defines.hpp
  dgt_dg.hpp
  dgt_field.hpp
  dgt_for_each.hpp
  dgt_grid3.hpp
  dgt_initialize.hpp
  dgt_lua.hpp
  dgt_lua_interface.hpp
  dgt_macros.hpp
  dgt_mesh.hpp
  dgt_partitioning.hpp
  dgt_print.hpp
  dgt_stl.hpp
  dgt_subgrid3.hpp
  dgt_tree.hpp
  dgt_vec.hpp
  dgt_vec3.hpp
  dgt_view.hpp
  dgt_when.hpp
)

set(DGT_SOURCES
  dgt_base64.cpp
  dgt_cartesian.cpp
  dgt_dg.cpp
  dgt_field.cpp
  dgt_initialize.cpp
  dgt_lua.cpp
  dgt_lua_interface.cpp
  dgt_mesh.cpp
  dgt_partitioning.cpp
  dgt_stl.cpp
  dgt_tree.cpp
  dgt_tree_vtu.cpp
  dgt_when.cpp
)

set_source_files_properties(${DGT_SOURCES} PROPERTIES LANGUAGE ${DGT_LANGUAGE})

add_library(dgtile ${DGT_SOURCES})
target_compile_features(dgtile PUBLIC cxx_std_17)
set_target_properties(dgtile PROPERTIES PUBLIC_HEADER "${DGT_HEADERS}")
set_target_properties(dgtile PROPERTIES ${DGT_LANGUAGE}_ARCHITECTURES "${DGT_ARCHITECTURES}")
target_include_directories(dgtile PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")
target_include_directories(dgtile PUBLIC "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
target_include_directories(dgtile PUBLIC "${LUA_INCLUDE_DIR}")
target_link_libraries(dgtile PUBLIC fmt::fmt)
target_link_libraries(dgtile PUBLIC Kokkos::kokkos)
target_link_libraries(dgtile PUBLIC ${LUA_LIBRARY})
target_link_libraries(dgtile PUBLIC mpicpp::mpicpp)
target_link_libraries(dgtile PUBLIC ZLIB::ZLIB)
install(TARGETS dgtile EXPORT dgtile-targets)
