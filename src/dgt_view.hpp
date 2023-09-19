#pragma once

#include <Kokkos_Core.hpp>

#include "dgt_macros.hpp"

namespace dgt {

template <class T>
using View = typename Kokkos::View<T, Kokkos::LayoutLeft>;

template <class T>
using HostView = typename Kokkos::View<T, Kokkos::LayoutLeft>::HostMirror;

#if defined(KOKKOS_ENABLE_CUDA)

template <class T>
using HostPinnedView = typename Kokkos::View<T,
      Kokkos::LayoutLeft,
      Kokkos::CudaHostPinnedSpace>;

template <class T>
using UnmanagedView = typename Kokkos::View<T,
      Kokkos::LayoutLeft,
      Kokkos::CudaHostPinnedSpace,
      Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

#elif defined(KOKKOS_ENABLE_HIP)

template <class T>
using HostPinnedView = typename Kokkos::View<T,
      Kokkos::LayoutLeft,
      Kokkos::HipHostPinnedSpace>;

template <class T>
using UnmanagedView = typename Kokkos::View<T,
      Kokkos::LayoutLeft,
      Kokkos::HipHostPinnedSpace,
      Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

#else

template <class T>
using HostPinnedView = typename Kokkos::View<T,
      Kokkos::LayoutLeft>;

template <class T>
using UnmanagedView = typename Kokkos::View<T,
      Kokkos::LayoutLeft,
      Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

#endif

template <class ViewT>
typename ViewT::value_type sum(ViewT v)
{
  using T = typename ViewT::value_type;
  T* data = v.data();
  T result = T(0);
  auto functor = [=] DGT_DEVICE (int const i, T& result) {
    result += data[i];
  };
  Kokkos::parallel_reduce("view_sum", v.size(), functor, result);
  return result;
}

}
