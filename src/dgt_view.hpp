#pragma once

#include <Kokkos_Core.hpp>

namespace dgt {

template <class T>
using View = typename Kokkos::View<T>;

template <class T>
using HostView = typename Kokkos::View<T>::HostMirror;

template <class T>
using UnmanagedView = typename Kokkos::View<T, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

#if defined(KOKKOS_ENABLE_CUDA)
template <class T>
using HostPinnedView = typename Kokkos::View<T, Kokkos::CudaHostPinnedSpace>;
#elif defined(KOKKOS_ENABLE_HIP)
template <class T>
using HostPinnedView = typename Kokkos::View<T, Kokkos::HipHostPinnedSpace>;
#else
template <class T>
using HostPinnedView = typename Kokkos::View<T>;
#endif

}
