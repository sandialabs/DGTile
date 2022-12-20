#pragma once

#include <stdexcept>

#include "Kokkos_Core.hpp"
#include "Kokkos_DualView.hpp"

#include "p3a_reduce.hpp"

namespace dgt {

template <class T>
using View = typename Kokkos::View<T, Kokkos::LayoutLeft>;

template <class T>
using HView = typename View<T>::HostMirror;

#if defined(KOKKOS_ENABLE_CUDA)
template <class T>
using HostPinnedView = typename Kokkos::View<T, Kokkos::LayoutLeft, Kokkos::CudaHostPinnedSpace>;
#elif defined(KOKKOS_ENABLE_HIP)
template <class T>
using HostPinnedView = typename Kokkos::View<T, Kokkos::LayoutLeft, Kokkos::HipHostPinnedSpace>;
#else
template <class T>
using HostPinnedView = HView<T>;
#endif

template <class InViewT, class OutViewT>
void resize(InViewT in, OutViewT& out) {
  int const n0 = in.extent(0);
  int const n1 = in.extent(1);
  int const n2 = in.extent(2);
  int const n3 = in.extent(3);
  int const n4 = in.extent(4);
  if (InViewT::Rank == 1) Kokkos::resize(out, n0);
  if (InViewT::Rank == 2) Kokkos::resize(out, n0, n1);
  if (InViewT::Rank == 3) Kokkos::resize(out, n0, n1, n2);
  if (InViewT::Rank == 4) Kokkos::resize(out, n0, n1, n2, n3);
  if (InViewT::Rank == 5) Kokkos::resize(out, n0, n1, n2, n3, n4);
}

template <class InViewT, class OutViewT>
void copy(InViewT in, OutViewT& out) {
  resize(in, out);
  Kokkos::deep_copy(out, in);
}

template <class ExecutionPolicy, class ViewT>
typename ViewT::value_type sum(ExecutionPolicy policy, ViewT v) {
  using T = typename ViewT::value_type;
  T* data = v.data();
  T constexpr identity_value = p3a::zero_value<T>();
  auto constexpr binary_op = p3a::adder<T>();
  auto functor = [=] P3A_HOST_DEVICE (int const i) {
    return data[i];
  };
  return p3a::transform_reduce(
      policy,
      p3a::counting_iterator(0),
      p3a::counting_iterator(int(v.size())),
      identity_value,
      binary_op,
      functor);
}

template <class ViewT>
void verify_extents(ViewT a, ViewT b) {
  for (int dim = 0; dim < ViewT::Rank; ++dim) {
    if (a.extent(dim) != b.extent(dim)) {
      throw std::runtime_error("incompatible views");
    }
  }
}

template <class ExecutionPolicy, class ViewT>
void axpby(
    ExecutionPolicy policy,
    ViewT r,
    typename ViewT::value_type const& a,
    ViewT x,
    typename ViewT::value_type const& b,
    ViewT y) {
  verify_extents(r, x);
  verify_extents(x, y);
  auto r_data = r.data();
  auto x_data = x.data();
  auto y_data = y.data();
  auto functor = [=] P3A_HOST_DEVICE (int const i) {
    r_data[i] = a*x_data[i] + b*y_data[i];
  };
  p3a::for_each(
      policy,
      p3a::counting_iterator(0),
      p3a::counting_iterator(int(r.size())),
      functor);
}

template <class ExecutionPolicy, class ViewT>
void xpaymz(
    ExecutionPolicy policy,
    ViewT r,
    ViewT x,
    typename ViewT::value_type const& a,
    ViewT y,
    ViewT z) {
  verify_extents(r, x);
  verify_extents(x, y);
  verify_extents(y, z);
  auto r_data = r.data();
  auto x_data = x.data();
  auto y_data = y.data();
  auto z_data = z.data();
  auto functor = [=] P3A_HOST_DEVICE (int const i) {
    r_data[i] = x_data[i] + a*(y_data[i] - z_data[i]);
  };
  p3a::for_each(
      policy,
      p3a::counting_iterator(0),
      p3a::counting_iterator(int(r.size())),
      functor);
}


}
