#pragma once

#include "p3a_static_vector.hpp"
#include "p3a_simd_view.hpp"

#include "dgt_basis.hpp"

namespace dgt {

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
device_simd<double> interp_scalar_intr(
    simd_view<double***> U, Basis const& b,
    int cell, int pt, int eq,
    device_simd_mask<double> const& mask) {
  device_simd<double> val = U.load(cell, eq, 0, mask) * b.phi_intr(pt, 0);
  for (int m = 1; m < b.nmodes; ++m) {
    val += U.load(cell, eq, m, mask) * b.phi_intr(pt, m);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
device_simd<double> interp_scalar_side(
    simd_view<double***> U, Basis const& b,
    int cell, int axis, int dir, int pt, int eq,
    device_simd_mask<double> const& mask) {
  device_simd<double> val = U.load(cell, eq, 0, mask) * b.phi_side(axis, dir, pt, 0);
  for (int m = 1; m < b.nmodes; ++m) {
    val += U.load(cell, eq, m, mask) * b.phi_side(axis, dir, pt, m);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
device_simd<double> interp_scalar_child_intr(
    simd_view<double***> U, Basis const& b,
    int cell, int which_child, int pt, int eq,
    device_simd_mask<double> const& mask) {
  device_simd<double> val = U.load(cell, eq, 0, mask) * b.phi_child_intr(which_child, pt, 0);
  for (int m = 1; m < b.nmodes; ++m) {
    val += U.load(cell, eq, m, mask) * b.phi_child_intr(which_child, pt, m);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
device_simd<double> interp_scalar_child_side(
    simd_view<double***> U, Basis const& b,
    int cell, int axis, int dir, int which_child, int pt, int eq,
    device_simd_mask<double> const& mask) {
  device_simd<double> val = U.load(cell, eq, 0, mask) * b.phi_child_side(axis, dir, which_child, pt, 0);
  for (int m = 1; m < b.nmodes; ++m) {
    val += U.load(cell, eq, m, mask) * b.phi_child_side(axis, dir, which_child, pt, m);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
device_simd<double> interp_scalar_fine(
    simd_view<double***> U, Basis const& b,
    int cell, int pt, int eq,
    device_simd_mask<double> const& mask) {
  device_simd<double> val = U.load(cell, eq, 0, mask) * b.phi_fine(pt, 0);
  for (int m = 1; m < b.nmodes; ++m) {
    val += U.load(cell, eq, m, mask) * b.phi_fine(pt, m);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
device_simd<double> interp_scalar_viz(
    simd_view<double***> U, Basis const& b,
    int cell, int pt, int eq,
    device_simd_mask<double> const& mask) {
  device_simd<double> val = U.load(cell, eq, 0, mask) * b.phi_viz(pt, 0);
  for (int m = 1; m < b.nmodes; ++m) {
    val += U.load(cell, eq, m, mask) * b.phi_viz(pt, m);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
device_simd<double> interp_scalar_eval(
    simd_view<double***> U, Basis const& b,
    int cell, int pt, int eq,
    device_simd_mask<double> const& mask) {
  device_simd<double> val = U.load(cell, eq, 0, mask) * b.phi_eval(pt, 0);
  for (int m = 1; m < b.nmodes; ++m) {
    val += U.load(cell, eq, m, mask) * b.phi_eval(pt, m);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<device_simd<double>> interp_vec3_intr(
    simd_view<double***> U, Basis const& b,
    int cell, int pt, int eq,
    device_simd_mask<double> const& mask) {
  vector3<device_simd<double>> val;
  for (int d = 0; d < DIMS; ++d) {
    val[d] = interp_scalar_intr(U, b, cell, pt, eq + d, mask);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<device_simd<double>> interp_vec3_side(
    simd_view<double***> U, Basis const& b,
    int cell, int axis, int dir, int pt, int eq,
    device_simd_mask<double> const& mask) {
  vector3<device_simd<double>> val;
  for (int d = 0; d < DIMS; ++d) {
    val[d] = interp_scalar_side(U, b, cell, axis, dir, pt, eq + d, mask);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<device_simd<double>> interp_vec3_child_intr(
    simd_view<double***> U, Basis const& b,
    int cell, int which_child, int pt, int eq,
    device_simd_mask<double> const& mask) {
  vector3<device_simd<double>> val;
  for (int d = 0; d < DIMS; ++d) {
    val[d] = interp_scalar_child_intr(U, b, cell, which_child, pt, eq + d, mask);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<device_simd<double>> interp_vec3_child_side(
    simd_view<double***> U, Basis const& b,
    int cell, int axis, int dir, int which_child, int pt, int eq,
    device_simd_mask<double> const& mask) {
  vector3<device_simd<double>> val;
  for (int d = 0; d < DIMS; ++d) {
    val[d] = interp_scalar_child_side(U, b, cell, axis, dir, which_child, pt, eq + d, mask);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<device_simd<double>> interp_vec3_fine(
    simd_view<double***> U, Basis const& b,
    int cell, int pt, int eq,
    device_simd_mask<double> const& mask) {
  vector3<device_simd<double>> val;
  for (int d = 0; d < DIMS; ++d) {
    val[d] = interp_scalar_fine(U, b, cell, pt, eq + d, mask);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<device_simd<double>> interp_vec3_viz(
    simd_view<double***> U, Basis const& b,
    int cell, int pt, int eq,
    device_simd_mask<double> const& mask) {
  vector3<device_simd<double>> val;
  for (int d = 0; d < DIMS; ++d) {
    val[d] = interp_scalar_viz(U, b, cell, pt, eq + d, mask);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<device_simd<double>> interp_vec3_eval(
    simd_view<double***> U, Basis const& b,
    int cell, int pt, int eq,
    device_simd_mask<double> const& mask) {
  vector3<device_simd<double>> val;
  for (int d = 0; d < DIMS; ++d) {
    val[d] = interp_scalar_eval(U, b, cell, pt, eq + d, mask);
  }
  return val;
}

template <int neq>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
static_vector<device_simd<double>, neq> interp_vec_intr(
    simd_view<double***> U, Basis const& b,
    int cell, int pt,
    device_simd_mask<double> const& mask) {
  static_vector<device_simd<double>, neq> val;
  for (int eq = 0; eq < neq; ++eq) {
    val[eq] = interp_scalar_intr(U, b, cell, pt, eq, mask);
  }
  return val;
}

template <int neq>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
static_vector<device_simd<double>, neq> interp_vec_side(
    simd_view<double***> U, Basis const& b,
    int cell, int axis, int dir, int pt,
    device_simd_mask<double> const& mask) {
  static_vector<device_simd<double>, neq> val;
  for (int eq = 0; eq < neq; ++eq) {
    val[eq] = interp_scalar_side(U, b, cell, axis, dir, pt, eq, mask);
  }
  return val;
}

template <int neq>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
static_vector<device_simd<double>, neq> interp_vec_child_intr(
    simd_view<double***> U, Basis const& b,
    int cell, int which_child, int pt,
    device_simd_mask<double> const& mask) {
  static_vector<device_simd<double>, neq> val;
  for (int eq = 0; eq < neq; ++eq) {
    val[eq] = interp_scalar_child_intr(U, b, cell, which_child, pt, eq, mask);
  }
  return val;
}

template <int neq>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
static_vector<device_simd<double>, neq> interp_vec_child_side(
    simd_view<double***> U, Basis const& b,
    int cell, int axis, int dir, int which_child, int pt,
    device_simd_mask<double> const& mask) {
  static_vector<device_simd<double>, neq> val;
  for (int eq = 0; eq < neq; ++eq) {
    val[eq] = interp_scalar_child_side(U, b, cell, axis, dir, which_child, pt, eq, mask);
  }
  return val;
}

template <int neq>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
static_vector<device_simd<double>, neq> interp_vec_fine(
    simd_view<double***> U, Basis const& b,
    int cell, int pt,
    device_simd_mask<double> const& mask) {
  static_vector<device_simd<double>, neq> val;
  for (int eq = 0; eq < neq; ++eq) {
    val[eq] = interp_scalar_fine(U, b, cell, pt, eq, mask);
  }
  return val;
}

template <int neq>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
static_vector<device_simd<double>, neq> interp_vec_viz(
    simd_view<double***> U, Basis const& b,
    int cell, int pt,
    device_simd_mask<double> const& mask) {
  static_vector<device_simd<double>, neq> val;
  for (int eq = 0; eq < neq; ++eq) {
    val[eq] = interp_scalar_viz(U, b, cell, pt, eq, mask);
  }
  return val;
}

template <int neq>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
static_vector<device_simd<double>, neq> interp_vec_eval(
    simd_view<double***> U, Basis const& b,
    int cell, int pt,
    device_simd_mask<double> const& mask) {
  static_vector<device_simd<double>, neq> val;
  for (int eq = 0; eq < neq; ++eq) {
    val[eq] = interp_scalar_eval(U, b, cell, pt, eq, mask);
  }
  return val;
}

template <int neq>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
static_vector<device_simd<double>, neq> gather_avg(
    simd_view<double***> U, int cell, device_simd_mask<double> const& mask) {
  static_vector<device_simd<double>, neq> val;
  for (int eq = 0; eq < neq; ++eq) {
    val[eq] = U.load(cell, eq, 0, mask);
  }
  return val;
}

}
