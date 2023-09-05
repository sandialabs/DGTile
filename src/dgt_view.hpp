#pragma once

#include <Kokkos_Core.hpp>

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

}
