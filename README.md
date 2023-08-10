# <img src="logo.png" width="150"> DGTile
> A portably performant discontinuous Galerkin adaptive mesh library

## What is this?

DGTile is a lightweight C++17 adaptive mesh library is meant to support
explicit discontinuous Galerkin applications on high performance computing
machines. The mesh supports block-structured adaptivity (isotropic refinement
and coarsening), where each block represents a Cartesian grid in one, two and
three spatial dimensions. Additionally, support for modal discontinuous
Galerkin basis functions is provided. Distributed memory parallelism (MPI) is
achieved by partitioning blocks over MPI ranks, and shared memory parallelism
is achieved by parallelizing functors over each block.

## Dependencies

DGTile has required dependencies on:

  * [fmt](https://github.com/fmtlib/fmt) - string formatting
  * [kokkos](https://github.com/kokkos/kokkos) - parallel algorithms
  * [zlib](https://github.com/zlib-ng/zlib-ng) - data compression
  * [googletest](https://github.com/google/googletest) - testing
  * [spdlog](https://github.com/gabime/spdlog) - fast C++ logging

##

At Sandia, DGTile is SCR 2806.0
