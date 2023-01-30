#pragma once

#include "p3a_static_vector.hpp"

#ifdef __GNUC__
#ifndef __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#endif

#include "p3a_simd_view.hpp"

#include "dgt_basis.hpp"

namespace dgt {

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::device_simd<double> interp_scalar_intr(
    p3a::simd_view<double***> U, Basis const& b,
    int cell, int pt, int eq,
    p3a::device_simd_mask<double> const& mask) {
  p3a::device_simd<double> val = U.load(cell, eq, 0, mask) * b.phi_intr(pt, 0);
  for (int m = 1; m < b.nmodes; ++m) {
    val += U.load(cell, eq, m, mask) * b.phi_intr(pt, m);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::device_simd<double> interp_scalar_side(
    p3a::simd_view<double***> U, Basis const& b,
    int cell, int axis, int dir, int pt, int eq,
    p3a::device_simd_mask<double> const& mask) {
  p3a::device_simd<double> val = U.load(cell, eq, 0, mask) * b.phi_side(axis, dir, pt, 0);
  for (int m = 1; m < b.nmodes; ++m) {
    val += U.load(cell, eq, m, mask) * b.phi_side(axis, dir, pt, m);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::device_simd<double> interp_scalar_child_intr(
    p3a::simd_view<double***> U, Basis const& b,
    int cell, int which_child, int pt, int eq,
    p3a::device_simd_mask<double> const& mask) {
  p3a::device_simd<double> val = U.load(cell, eq, 0, mask) * b.phi_child_intr(which_child, pt, 0);
  for (int m = 1; m < b.nmodes; ++m) {
    val += U.load(cell, eq, m, mask) * b.phi_child_intr(which_child, pt, m);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::device_simd<double> interp_scalar_child_side(
    p3a::simd_view<double***> U, Basis const& b,
    int cell, int axis, int dir, int which_child, int pt, int eq,
    p3a::device_simd_mask<double> const& mask) {
  p3a::device_simd<double> val = U.load(cell, eq, 0, mask) * b.phi_child_side(axis, dir, which_child, pt, 0);
  for (int m = 1; m < b.nmodes; ++m) {
    val += U.load(cell, eq, m, mask) * b.phi_child_side(axis, dir, which_child, pt, m);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::device_simd<double> interp_scalar_fine(
    p3a::simd_view<double***> U, Basis const& b,
    int cell, int pt, int eq,
    p3a::device_simd_mask<double> const& mask) {
  p3a::device_simd<double> val = U.load(cell, eq, 0, mask) * b.phi_fine(pt, 0);
  for (int m = 1; m < b.nmodes; ++m) {
    val += U.load(cell, eq, m, mask) * b.phi_fine(pt, m);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::device_simd<double> interp_scalar_eval(
    p3a::simd_view<double***> U, Basis const& b,
    int cell, int pt, int eq,
    p3a::device_simd_mask<double> const& mask) {
  p3a::device_simd<double> val = U.load(cell, eq, 0, mask) * b.phi_eval(pt, 0);
  for (int m = 1; m < b.nmodes; ++m) {
    val += U.load(cell, eq, m, mask) * b.phi_eval(pt, m);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::device_simd<double> interp_scalar_corner(
    p3a::simd_view<double***> U, Basis const& b,
    int cell, int pt, int eq,
    p3a::device_simd_mask<double> const& mask) {
  p3a::device_simd<double> val = U.load(cell, eq, 0, mask) * b.phi_corner(pt, 0);
  for (int m = 1; m < b.nmodes; ++m) {
    val += U.load(cell, eq, m, mask) * b.phi_corner(pt, m);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::vector3<p3a::device_simd<double>> interp_vec3_intr(
    p3a::simd_view<double***> U, Basis const& b,
    int cell, int pt, int eq,
    p3a::device_simd_mask<double> const& mask) {
  p3a::vector3<p3a::device_simd<double>> val;
  for (int d = 0; d < DIMS; ++d) {
    val[d] = interp_scalar_intr(U, b, cell, pt, eq + d, mask);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::vector3<p3a::device_simd<double>> interp_vec3_side(
    p3a::simd_view<double***> U, Basis const& b,
    int cell, int axis, int dir, int pt, int eq,
    p3a::device_simd_mask<double> const& mask) {
  p3a::vector3<p3a::device_simd<double>> val;
  for (int d = 0; d < DIMS; ++d) {
    val[d] = interp_scalar_side(U, b, cell, axis, dir, pt, eq + d, mask);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::vector3<p3a::device_simd<double>> interp_vec3_child_intr(
    p3a::simd_view<double***> U, Basis const& b,
    int cell, int which_child, int pt, int eq,
    p3a::device_simd_mask<double> const& mask) {
  p3a::vector3<p3a::device_simd<double>> val;
  for (int d = 0; d < DIMS; ++d) {
    val[d] = interp_scalar_child_intr(U, b, cell, which_child, pt, eq + d, mask);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::vector3<p3a::device_simd<double>> interp_vec3_child_side(
    p3a::simd_view<double***> U, Basis const& b,
    int cell, int axis, int dir, int which_child, int pt, int eq,
    p3a::device_simd_mask<double> const& mask) {
  p3a::vector3<p3a::device_simd<double>> val;
  for (int d = 0; d < DIMS; ++d) {
    val[d] = interp_scalar_child_side(U, b, cell, axis, dir, which_child, pt, eq + d, mask);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::vector3<p3a::device_simd<double>> interp_vec3_fine(
    p3a::simd_view<double***> U, Basis const& b,
    int cell, int pt, int eq,
    p3a::device_simd_mask<double> const& mask) {
  p3a::vector3<p3a::device_simd<double>> val;
  for (int d = 0; d < DIMS; ++d) {
    val[d] = interp_scalar_fine(U, b, cell, pt, eq + d, mask);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::vector3<p3a::device_simd<double>> interp_vec3_eval(
    p3a::simd_view<double***> U, Basis const& b,
    int cell, int pt, int eq,
    p3a::device_simd_mask<double> const& mask) {
  p3a::vector3<p3a::device_simd<double>> val;
  for (int d = 0; d < DIMS; ++d) {
    val[d] = interp_scalar_eval(U, b, cell, pt, eq + d, mask);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::vector3<p3a::device_simd<double>> interp_vec3_corner(
    p3a::simd_view<double***> U, Basis const& b,
    int cell, int pt, int eq,
    p3a::device_simd_mask<double> const& mask) {
  p3a::vector3<p3a::device_simd<double>> val;
  for (int d = 0; d < DIMS; ++d) {
    val[d] = interp_scalar_corner(U, b, cell, pt, eq + d, mask);
  }
  return val;
}

template <int neq>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::static_vector<p3a::device_simd<double>, neq> interp_vec_intr(
    p3a::simd_view<double***> U, Basis const& b,
    int cell, int pt,
    p3a::device_simd_mask<double> const& mask) {
  p3a::static_vector<p3a::device_simd<double>, neq> val;
  for (int eq = 0; eq < neq; ++eq) {
    val[eq] = interp_scalar_intr(U, b, cell, pt, eq, mask);
  }
  return val;
}

template <int neq>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::static_vector<p3a::device_simd<double>, neq> interp_vec_side(
    p3a::simd_view<double***> U, Basis const& b,
    int cell, int axis, int dir, int pt,
    p3a::device_simd_mask<double> const& mask) {
  p3a::static_vector<p3a::device_simd<double>, neq> val;
  for (int eq = 0; eq < neq; ++eq) {
    val[eq] = interp_scalar_side(U, b, cell, axis, dir, pt, eq, mask);
  }
  return val;
}

template <int neq>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::static_vector<p3a::device_simd<double>, neq> interp_vec_child_intr(
    p3a::simd_view<double***> U, Basis const& b,
    int cell, int which_child, int pt,
    p3a::device_simd_mask<double> const& mask) {
  p3a::static_vector<p3a::device_simd<double>, neq> val;
  for (int eq = 0; eq < neq; ++eq) {
    val[eq] = interp_scalar_child_intr(U, b, cell, which_child, pt, eq, mask);
  }
  return val;
}

template <int neq>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::static_vector<p3a::device_simd<double>, neq> interp_vec_child_side(
    p3a::simd_view<double***> U, Basis const& b,
    int cell, int axis, int dir, int which_child, int pt,
    p3a::device_simd_mask<double> const& mask) {
  p3a::static_vector<p3a::device_simd<double>, neq> val;
  for (int eq = 0; eq < neq; ++eq) {
    val[eq] = interp_scalar_child_side(U, b, cell, axis, dir, which_child, pt, eq, mask);
  }
  return val;
}

template <int neq>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::static_vector<p3a::device_simd<double>, neq> interp_vec_fine(
    p3a::simd_view<double***> U, Basis const& b,
    int cell, int pt,
    p3a::device_simd_mask<double> const& mask) {
  p3a::static_vector<p3a::device_simd<double>, neq> val;
  for (int eq = 0; eq < neq; ++eq) {
    val[eq] = interp_scalar_fine(U, b, cell, pt, eq, mask);
  }
  return val;
}

template <int neq>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::static_vector<p3a::device_simd<double>, neq> interp_vec_eval(
    p3a::simd_view<double***> U, Basis const& b,
    int cell, int pt,
    p3a::device_simd_mask<double> const& mask) {
  p3a::static_vector<p3a::device_simd<double>, neq> val;
  for (int eq = 0; eq < neq; ++eq) {
    val[eq] = interp_scalar_eval(U, b, cell, pt, eq, mask);
  }
  return val;
}

template <int neq>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::static_vector<p3a::device_simd<double>, neq> interp_vec_corner(
    p3a::simd_view<double***> U, Basis const& b,
    int cell, int pt,
    p3a::device_simd_mask<double> const& mask) {
  p3a::static_vector<p3a::device_simd<double>, neq> val;
  for (int eq = 0; eq < neq; ++eq) {
    val[eq] = interp_scalar_corner(U, b, cell, pt, eq, mask);
  }
  return val;
}

template <int neq>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
p3a::static_vector<p3a::device_simd<double>, neq> gather_avg(
    p3a::simd_view<double***> U, int cell, p3a::device_simd_mask<double> const& mask) {
  p3a::static_vector<p3a::device_simd<double>, neq> val;
  int const m = 0;
  for (int eq = 0; eq < neq; ++eq) {
    val[eq] = U.load(cell, eq, m, mask);
  }
  return val;
}

}

#ifdef __GNUC__
#ifndef __clang__
#pragma GCC diagnostic pop
#endif
#endif
