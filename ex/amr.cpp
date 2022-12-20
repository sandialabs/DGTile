#include <stdexcept>

#include "caliper/cali.h"

#include "dgt_amr.hpp"
#include "dgt_grid.hpp"
#include "dgt_marks.hpp"
#include "dgt_mesh.hpp"

#include "hydro.hpp"

namespace hydro {

static void refine_first_block_initial(Mesh& mesh) {
  CALI_CXX_MARK_FUNCTION;
  refine(mesh.dim(), mesh.leaves()[0]);
  mesh.rebuild();
}

static void refine_rt_initial(Mesh& mesh) {
  CALI_CXX_MARK_FUNCTION;
  refine(mesh.dim(), mesh.leaves()[1]);
  mesh.rebuild();
  for (Node* leaf : mesh.leaves()) {
    refine(mesh.dim(), leaf);
  }
  mesh.rebuild();
  refine(mesh.dim(), mesh.leaves()[6]);
  refine(mesh.dim(), mesh.leaves()[7]);
  refine(mesh.dim(), mesh.leaves()[10]);
  refine(mesh.dim(), mesh.leaves()[11]);
  refine(mesh.dim(), mesh.leaves()[12]);
  refine(mesh.dim(), mesh.leaves()[13]);
  refine(mesh.dim(), mesh.leaves()[16]);
  refine(mesh.dim(), mesh.leaves()[17]);
  mesh.rebuild();
}

static void do_clockwise_amr(State& state) {
  CALI_CXX_MARK_FUNCTION;
  static int ctr = 1;
  std::vector<int8_t> marks(state.mesh.leaves().size(), REMAIN);
  double const freq = state.in.amr_frequency;
  if (state.t >= freq * ctr) {
    if (state.mesh.comm()->rank() == 0) {
      std::cout << "adapting\n";
    }
    if (ctr == 1) {
      marks[0] = COARSEN;
      marks[1] = COARSEN;
      marks[2] = COARSEN;
      marks[3] = COARSEN;
      marks[4] = REFINE;
      modify(state.mesh, marks);
    }
    if (ctr == 2) {
      marks[1] = COARSEN;
      marks[2] = COARSEN;
      marks[3] = COARSEN;
      marks[4] = COARSEN;
      marks[6] = REFINE;
      modify(state.mesh, marks);
    }
    if (ctr == 3) {
      marks[3] = COARSEN;
      marks[4] = COARSEN;
      marks[5] = COARSEN;
      marks[6] = COARSEN;
      marks[2] = REFINE;
      modify(state.mesh, marks);
    }
    if (ctr == 4) {
      marks[2] = COARSEN;
      marks[3] = COARSEN;
      marks[4] = COARSEN;
      marks[5] = COARSEN;
      marks[0] = REFINE;
      modify(state.mesh, marks);
    }
    ctr++;
  }
}

static void do_rt_amr(State& state) {
  CALI_CXX_MARK_FUNCTION;
  double const tol = 1.e-3;
  static int ctr = 1;
  double const freq = state.in.amr_frequency;
  if (state.t >= freq*ctr) {
    if (state.mesh.comm()->rank() == 0) {
      std::cout << "adapting\n";
    }
    Mesh& mesh = state.mesh;
    int const dim = mesh.dim();
    std::vector<int8_t> marks(mesh.leaves().size(), REMAIN);
    for (Node* leaf : mesh.owned_leaves()) {
      if (leaf->pt().depth == 6) continue;
      device_array<int8_t> mark_vector;
      mark_vector.resize(1, REMAIN);
      Block& block = leaf->block;
      grid3 const cell_grid = generalize(block.cell_grid());
      double const detJ = block.cell_detJ();
      double const volume = get_volume(dim, block.domain().extents());
      View<double***> U = block.soln(0);
      auto f = [=] P3A_DEVICE (vector3<int> const& cell_ijk) {
        int const cell = cell_grid.index(cell_ijk);
        double val = 0.;
        for (int axis = 0; axis < dim; ++axis) {
          val += std::abs(U(cell, RH, 1 + axis)) * detJ;
        }
        return val;
      };
      double constexpr identity_value = zero_value<double>();
      auto constexpr binary_op = adder<double>();
      double const result = transform_reduce(
          execution::par, cell_grid, identity_value, binary_op, f);
      double const val = result / volume;
      if (val > tol) {
        marks[block.id()] = REFINE;
      }
    }
    modify(state.mesh, marks);
    ctr++;
  }
}

static void do_debug_amr(State& state) {
  CALI_CXX_MARK_FUNCTION;
  static int ctr = 0;
  std::vector<int8_t> marks(state.mesh.leaves().size(), REMAIN);
  if (ctr == 0) {
    if (state.mesh.comm()->rank() == 0) {
      std::cout << "adapting\n";
    }
    marks[3] = REFINE;
    modify(state.mesh, marks);
    ctr++;
  }
}

void do_initial_amr(State& state) {
  CALI_CXX_MARK_FUNCTION;
  std::string const amr  = state.in.init_amr;
  if (amr == "") {
    return;
  } else if (amr == "first_block") {
    refine_first_block_initial(state.mesh);
  } else if (amr == "rt") {
    refine_rt_initial(state.mesh);
  } else {
    throw std::runtime_error("invalid init_amr");
  }
}

void do_amr(State& state) {
  CALI_CXX_MARK_FUNCTION;
  std::string const amr = state.in.amr;
  if (amr == "") {
    return;
  } else if (amr == "clockwise") {
    do_clockwise_amr(state);
  } else if (amr == "rt") {
    do_rt_amr(state);
  } else if (amr == "debug") {
    do_debug_amr(state);
  } else {
    throw std::runtime_error("invalid amr");
  }
}

}
