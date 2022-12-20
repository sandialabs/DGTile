# DGTile

Portably performant discontinuous Galerkin adaptive mesh library

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

  * [p3a](https://github.com/sandialabs/p3a) - structured grid support
    * [kokkos](https://github.com/kokkos/kokkos) - parallel algorithms
    * [mpicpp](https://github.com/sandialabs/mpicpp) - MPI interface
  * [caliper](https://github.com/LLNL/Caliper) - performance profiling
  * [zlib](https://github.com/zlib-ng/zlib-ng) - data compression

DGTile has an optional dependency on:

  * [googletest](https://github.com/google/googletest)

if `-DBUILD_TESTING=ON`.

##

At Sandia, DGTile is SCR 2806.0
