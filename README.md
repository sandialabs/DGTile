<img src="logo.png" width="200">

# DGTile
> A portably performant discontinuous Galerkin adaptive mesh library

## What is this?

DGTile is a C++17 adaptive mesh library meant to support explicit
discontinuous Galerkin finite element applications on high performance
computers.

## Dependencies

DGTile has required dependencies on:

  * [fmt](https://github.com/fmtlib/fmt) - string formatting
  * [kokkos](https://github.com/kokkos/kokkos) - shared memory parallelism
  * [googletest](https://github.com/google/googletest) - testing
  * [lua](https://gitlab.com/codelibre/lua/lua-cmake) - interfacing
  * [mpi](https://www.open-mpi.org/) - distributed memory parallelism
  * [mpicpp](https://github.com/sandialabs/mpicpp) - MPI C++ interfacing
  * [zlib](https://github.com/zlib-ng/zlib-ng) - data compression

## Documentation

See the [documentation](https://sandialabs.github.io/DGTile/index.html)
for information about:

  * deployment (compiling and installing the code)
  * mathematical, numerical, and computation theory
  * usage (how the code can be used in applications)
  * code performance

##

At Sandia, DGTile is SCR 2806.0
