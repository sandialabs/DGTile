cmake_minimum_required(VERSION 3.12.0...3.18.0)

project(DGTile VERSION 0.1.0 LANGUAGES CXX C)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

find_package(p3a REQUIRED)
find_package(ZLIB REQUIRED)
find_package(caliper REQUIRED)

option(dgtile_WORKAROUND_GCC_FILESYSTEM "link to std++fs" OFF)

if (BUILD_TESTING)
  enable_testing()
  find_package(GTest REQUIRED)
endif()

if (BUILD_EXAMPLE)
  enable_testing()
endif()

enable_language(${p3a_LANGUAGE})
if (p3a_ENABLE_CUDA)
  find_package(CUDAToolkit REQUIRED)
endif()

set(DGT_HEADERS
  dgt_amr.hpp
  dgt_basis.hpp
  dgt_block.hpp
  dgt_border.hpp
  dgt_defines.hpp
  dgt_field.hpp
  dgt_file.hpp
  dgt_grid.hpp
  dgt_interp.hpp
  dgt_interp_simd.hpp
  dgt_library.hpp
  dgt_marks.hpp
  dgt_mesh.hpp
  dgt_message.hpp
  dgt_point.hpp
  dgt_print.hpp
  dgt_spatial.hpp
  dgt_tree.hpp
  dgt_views.hpp
)

set(DGT_SOURCES
  dgt_amr.cpp
  dgt_ascii.cpp
  dgt_base64.cpp
  dgt_basis.cpp
  dgt_block.cpp
  dgt_binary.cpp
  dgt_border.cpp
  dgt_field.cpp
  dgt_file.cpp
  dgt_grid.cpp
  dgt_library.cpp
  dgt_marks.cpp
  dgt_mesh.cpp
  dgt_tree.cpp
  dgt_vtk.cpp
)

set_source_files_properties(${DGT_SOURCES} PROPERTIES LANGUAGE ${p3a_LANGUAGE})

add_library(dgtile ${DGT_SOURCES})

target_include_directories(
  dgtile
  PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")

target_compile_features(dgtile PUBLIC cxx_std_17)

set_target_properties(dgtile PROPERTIES PUBLIC_HEADER "${DGT_HEADERS}")
set_target_properties(dgtile PROPERTIES ${p3a_LANGUAGE}_ARCHITECTURES "${p3a_ARCHITECTURES}")

target_link_libraries(dgtile PUBLIC MPI::MPI_CXX)
target_link_libraries(dgtile PUBLIC mpicpp::mpicpp)
target_link_libraries(dgtile PUBLIC p3a::p3a)
target_link_libraries(dgtile PUBLIC ZLIB::ZLIB)
target_link_libraries(dgtile PUBLIC caliper)

if (dgtile_WORKAROUND_GCC_FILESYSTEM)
  target_link_libraries(dgtile PUBLIC stdc++fs)
endif()

install(
  TARGETS dgtile
  EXPORT dgtile-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if (BUILD_TESTING)
  add_subdirectory(unit)
endif()

if (BUILD_EXAMPLE)
  add_subdirectory(ex)
endif()

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/dgtile-config.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/dgtile
)

write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/dgtile-config-version.cmake"
  COMPATIBILITY AnyNewerVersion)

install(FILES
  "${PROJECT_BINARY_DIR}/dgtile-config.cmake"
  "${PROJECT_BINARY_DIR}/dgtile-config-version.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/dgtile)

install(
  EXPORT dgtile-targets
  NAMESPACE dgtile::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/dgtile)
