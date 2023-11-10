cmake_minimum_required(VERSION 3.12.0...3.18.0)

project(DGTile VERSION 0.1.0 LANGUAGES CXX C)

enable_testing()
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(DGT_ENABLE_CONDUIT "Use conduit for file IO" OFF)

set(DGT_ARCHITECTURES "" CACHE STRING "GPU architecture")

find_package(fmt REQUIRED)
find_package(GTest REQUIRED)
find_package(Kokkos REQUIRED COMPONENTS separable_compilation)
find_package(Lua REQUIRED)
find_package(MPI REQUIRED)
find_package(mpicpp REQUIRED)
find_package(ZLIB REQUIRED)

if(DGT_ENABLE_CONDUIT)
  find_package(Conduit REQUIRED)
endif()

set(DGT_LANGUAGE CXX)
if(Kokkos_COMPILE_LANGUAGE)
  set(DGT_LANGUAGE "${Kokkos_COMPILE_LANGUAGE}")
endif()
enable_language(${DGT_LANGUAGE})

add_subdirectory(src)
add_subdirectory(test)

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
