#pragma once

#include "p3a_static_vector.hpp"

#include "dgt_basis.hpp"

namespace dgt {

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
double interp_scalar_intr(
    View<double***> U, Basis const& b,
    int cell, int pt, int eq) {
  double val = U(cell, eq, 0) * b.phi_intr(pt, 0);
  for (int m = 1; m < b.nmodes; ++m) {
    val += U(cell, eq, m) * b.phi_intr(pt, m);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
double interp_scalar_side(
    View<double***> U, Basis const& b,
    int cell, int axis, int dir, int pt, int eq) {
  double val = U(cell, eq, 0) * b.phi_side(axis, dir, pt, 0);
  for (int m = 1; m < b.nmodes; ++m) {
    val += U(cell, eq, m) * b.phi_side(axis, dir, pt, m);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
double interp_scalar_child_intr(
    View<double***> U, Basis const& b,
    int cell, int which_child, int pt, int eq) {
  double val = U(cell, eq, 0) * b.phi_child_intr(which_child, pt, 0);
  for (int m = 1; m < b.nmodes; ++m) {
    val += U(cell, eq, m) * b.phi_child_intr(which_child, pt, m);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
double interp_scalar_child_side(
    View<double***> U, Basis const& b,
    int cell, int axis, int dir, int which_child, int pt, int eq) {
  double val = U(cell, eq, 0) * b.phi_child_side(axis, dir, which_child, pt, 0);
  for (int m = 1; m < b.nmodes; ++m) {
    val += U(cell, eq, m) * b.phi_child_side(axis, dir, which_child, pt, m);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
double interp_scalar_fine(
    View<double***> U, Basis const& b,
    int cell, int pt, int eq) {
  double val = U(cell, eq, 0) * b.phi_fine(pt, 0);
  for (int m = 1; m < b.nmodes; ++m) {
    val += U(cell, eq, m) * b.phi_fine(pt, m);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
double interp_scalar_viz(
    View<double***> U, Basis const& b,
    int cell, int pt, int eq) {
  double val = U(cell, eq, 0) * b.phi_viz(pt, 0);
  for (int m = 1; m < b.nmodes; ++m) {
    val += U(cell, eq, m) * b.phi_viz(pt, m);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
double interp_scalar_eval(
    View<double***> U, Basis const& b,
    int cell, int pt, int eq) {
  double val = U(cell, eq, 0) * b.phi_eval(pt, 0);
  for (int m = 1; m < b.nmodes; ++m) {
    val += U(cell, eq, m) * b.phi_eval(pt, m);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<double> interp_vec3_intr(
    View<double***> U, Basis const& b,
    int cell, int pt, int eq) {
  vector3<double> val;
  for (int d = 0; d < DIMS; ++d) {
    val[d] = interp_scalar_intr(U, b, cell, pt, eq + d);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<double> interp_vec3_side(
    View<double***> U, Basis const& b,
    int cell, int axis, int dir, int pt, int eq) {
  vector3<double> val;
  for (int d = 0; d < DIMS; ++d) {
    val[d] = interp_scalar_side(U, b, cell, axis, dir, pt, eq + d);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<double> interp_vec3_child_intr(
    View<double***> U, Basis const& b,
    int cell, int which_child, int pt, int eq) {
  vector3<double> val;
  for (int d = 0; d < DIMS; ++d) {
    val[d] = interp_scalar_child_intr(U, b, cell, which_child, pt, eq + d);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<double> interp_vec3_child_side(
    View<double***> U, Basis const& b,
    int cell, int axis, int dir, int which_child, int pt, int eq) {
  vector3<double> val;
  for (int d = 0; d < DIMS; ++d) {
    val[d] = interp_scalar_child_side(U, b, cell, axis, dir, which_child, pt, eq + d);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<double> interp_vec3_fine(
    View<double***> U, Basis const& b,
    int cell, int pt, int eq) {
  vector3<double> val;
  for (int d = 0; d < DIMS; ++d) {
    val[d] = interp_scalar_fine(U, b, cell, pt, eq + d);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<double> interp_vec3_viz(
    View<double***> U, Basis const& b,
    int cell, int pt, int eq) {
  vector3<double> val;
  for (int d = 0; d < DIMS; ++d) {
    val[d] = interp_scalar_viz(U, b, cell, pt, eq + d);
  }
  return val;
}

[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
vector3<double> interp_vec3_eval(
    View<double***> U, Basis const& b,
    int cell, int pt, int eq) {
  vector3<double> val;
  for (int d = 0; d < DIMS; ++d) {
    val[d] = interp_scalar_eval(U, b, cell, pt, eq + d);
  }
  return val;
}

template <int neq>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
static_vector<double, neq> interp_vec_intr(
    View<double***> U, Basis const& b,
    int cell, int pt) {
  static_vector<double, neq> val;
  for (int eq = 0; eq < neq; ++eq) {
    val[eq] = interp_scalar_intr(U, b, cell, pt, eq);
  }
  return val;
}

template <int neq>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
static_vector<double, neq> interp_vec_side(
    View<double***> U, Basis const& b,
    int cell, int axis, int dir, int pt) {
  static_vector<double, neq> val;
  for (int eq = 0; eq < neq; ++eq) {
    val[eq] = interp_scalar_side(U, b, cell, axis, dir, pt, eq);
  }
  return val;
}

template <int neq>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
static_vector<double, neq> interp_vec_child_intr(
    View<double***> U, Basis const& b,
    int cell, int which_child, int pt) {
  static_vector<double, neq> val;
  for (int eq = 0; eq < neq; ++eq) {
    val[eq] = interp_scalar_child_intr(U, b, cell, which_child, pt, eq);
  }
  return val;
}

template <int neq>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
static_vector<double, neq> interp_vec_child_side(
    View<double***> U, Basis const& b,
    int cell, int axis, int dir, int which_child, int pt) {
  static_vector<double, neq> val;
  for (int eq = 0; eq < neq; ++eq) {
    val[eq] = interp_scalar_child_side(U, b, cell, axis, dir, which_child, pt, eq);
  }
  return val;
}

template <int neq>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
static_vector<double, neq> interp_vec_fine(
    View<double***> U, Basis const& b,
    int cell, int pt) {
  static_vector<double, neq> val;
  for (int eq = 0; eq < neq; ++eq) {
    val[eq] = interp_scalar_fine(U, b, cell, pt, eq);
  }
}

template <int neq>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
static_vector<double, neq> interp_vec_viz(
    View<double***> U, Basis const& b,
    int cell, int pt) {
  static_vector<double, neq> val;
  for (int eq = 0; eq < neq; ++eq) {
    val[eq] = interp_scalar_viz(U, b, cell, pt, eq);
  }
  return val;
}

template <int neq>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
static_vector<double, neq> interp_vec_eval(
    View<double***> U, Basis const& b,
    int cell, int pt) {
  static_vector<double, neq> val;
  for (int eq = 0; eq < neq; ++eq) {
    val[eq] = interp_scalar_eval(U, b, cell, pt, eq);
  }
  return val;
}

template <int neq>
[[nodiscard]] P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline
static_vector<double, neq> gather_avg(View<double***> U, int cell) {
  static_vector<double, neq> val;
  for (int eq = 0; eq < neq; ++eq) {
    val[eq] = U(cell, eq, 0);
  }
  return val;
}

}
