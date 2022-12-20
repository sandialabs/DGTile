#pragma once

#include <stdexcept>

#include "Kokkos_Core.hpp"
#include "Kokkos_DualView.hpp"

namespace dgt {

template <class T>
using View = typename Kokkos::View<T, Kokkos::LayoutLeft>;

template <class A, class B>
void copy(A in, B& out) {
  int const n0 = in.extent(0);
  int const n1 = in.extent(1);
  int const n2 = in.extent(2);
  int const n3 = in.extent(3);
  int const n4 = in.extent(4);
  if (A::Rank == 1) Kokkos::resize(out, n0);
  if (A::Rank == 2) Kokkos::resize(out, n0, n1);
  if (A::Rank == 3) Kokkos::resize(out, n0, n1, n2);
  if (A::Rank == 4) Kokkos::resize(out, n0, n1, n2, n3);
  if (A::Rank == 5) Kokkos::resize(out, n0, n1, n2, n3, n4);
  Kokkos::deep_copy(out, in);
}

template <class ViewT>
typename ViewT::value_type sum(ViewT view) {
  using T = typename ViewT::value_type;
  T val = T(0);
  T* data = view.data();
  auto f = KOKKOS_LAMBDA (int const i, T& val) {
    val += data[i];
  };
  int const size = view.size();
  Kokkos::parallel_reduce(size, f, Kokkos::Sum<T>(val));
  Kokkos::fence();
  return val;
}

template <class ViewT>
void verify_extents(ViewT a, ViewT b) {
  for (int dim = 0; dim < ViewT::Rank; ++dim) {
    if (a.extent(dim) != b.extent(dim)) {
      throw std::runtime_error("incompatible views");
    }
  }
}

template <class ViewT>
void axpby(
    ViewT r,
    typename ViewT::value_type const& a,
    ViewT x,
    typename ViewT::value_type const& b,
    ViewT y) {
  verify_extents(r,x);
  verify_extents(x,y);
  typename ViewT::value_type* rdata = r.data();
  typename ViewT::value_type* xdata = x.data();
  typename ViewT::value_type* ydata = y.data();
  auto f = KOKKOS_LAMBDA (int const i) {
    rdata[i] = a*xdata[i] + b*ydata[i];
  };
  Kokkos::parallel_for("axby", r.size(), f);
  Kokkos::fence();
}

}
