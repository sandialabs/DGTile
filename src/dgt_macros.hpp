#pragma once

#include <Kokkos_Macros.hpp>

#define DGT_ALWAYS_INLINE __attribute__((always_inline))
#define DGT_NEVER_INLINE __attribute__((noinline))

#ifdef KOKKOS_ENABLE_CUDA
#ifndef __CUDACC_EXTENDED_LAMBDA__
#error "please recompile with -expt-extended-lambda"
#endif
#endif

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
#define DGT_HOST __host__
#define DGT_DEVICE __device__
#else
#define DGT_HOST
#define DGT_DEVICE
#endif

#define DGT_HOST_DEVICE DGT_HOST DGT_DEVICE

#define DGT_METHOD [[nodiscard]] DGT_ALWAYS_INLINE DGT_HOST DGT_DEVICE
#define DGT_VOID_METHOD DGT_ALWAYS_INLINE DGT_HOST DGT_DEVICE
