#include <stdexcept>

#include "caliper/cali.h"

#include "dgt_amr.hpp"
#include "dgt_grid.hpp"
#include "dgt_marks.hpp"
#include "dgt_mesh.hpp"

#include "hydro.hpp"

namespace hydro {

static void refine_blocks(Mesh& mesh, std::vector<int> const& blocks) {
  for (int b : blocks) {
    refine(mesh.dim(), mesh.leaves()[b]);
  }
}

static void refine_first_block_initial(Mesh& mesh) {
  CALI_CXX_MARK_FUNCTION;
  refine(mesh.dim(), mesh.leaves()[0]);
  mesh.rebuild();
}

static void refine_sod_initial(Mesh& mesh) {
  CALI_CXX_MARK_FUNCTION;
  refine_blocks(mesh, {1,2});
  mesh.rebuild();
  refine_blocks(mesh, {2,4,5,7});
  mesh.rebuild();
  refine_blocks(mesh, {3,5,8,10,11,13,16,18});
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
  refine_blocks(mesh, {6,7,10,11,12,13,16,17});
  mesh.rebuild();
}

static void do_clockwise_amr(State& state) {
  CALI_CXX_MARK_FUNCTION;
  static int ctr = 1;
  std::vector<int8_t> marks(state.mesh.leaves().size(), dgt::REMAIN);
  double const freq = state.in.amr_frequency;
  if (state.t >= freq * ctr) {
    if (state.mesh.comm()->rank() == 0) {
      std::cout << "adapting\n";
    }
    if (ctr == 1) {
      marks[0] = dgt::COARSEN;
      marks[1] = dgt::COARSEN;
      marks[2] = dgt::COARSEN;
      marks[3] = dgt::COARSEN;
      marks[4] = dgt::REFINE;
      modify(state.mesh, marks);
    }
    if (ctr == 2) {
      marks[1] = dgt::COARSEN;
      marks[2] = dgt::COARSEN;
      marks[3] = dgt::COARSEN;
      marks[4] = dgt::COARSEN;
      marks[6] = dgt::REFINE;
      modify(state.mesh, marks);
    }
    if (ctr == 3) {
      marks[3] = dgt::COARSEN;
      marks[4] = dgt::COARSEN;
      marks[5] = dgt::COARSEN;
      marks[6] = dgt::COARSEN;
      marks[2] = dgt::REFINE;
      modify(state.mesh, marks);
    }
    if (ctr == 4) {
      marks[2] = dgt::COARSEN;
      marks[3] = dgt::COARSEN;
      marks[4] = dgt::COARSEN;
      marks[5] = dgt::COARSEN;
      marks[0] = dgt::REFINE;
      modify(state.mesh, marks);
    }
    ctr++;
  }
}

#if 0
static void do_sod_amr(State& state) {
  CALI_CXX_MARK_FUNCTION;
  static int ctr = 1;
  double const freq = state.in.amr_frequency;
  if (state.t >= freq*ctr) {
    if (state.mesh.com()->rank() == 0) {
      std::cout << "adapting\n";
    }
    Mesh& mesh = state.mesh;
    int const dim = mesh.dim();
    std::vector<int8_t> marks(mesh.leaves().size(), dgt::REMAIN);
    if (leaf->pt().deth == 5) continue;
    Block& block = leaff->block;
    p3a::grid3 const cell_grid = dgt::generalize(block.cell_grid());
    double const detJ = block.cell_detJ();
    ctr++;
  }
}
#endif

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
    std::vector<int8_t> marks(mesh.leaves().size(), dgt::REMAIN);
    for (Node* leaf : mesh.owned_leaves()) {
      if (leaf->pt().depth == 6) continue;
      Block& block = leaf->block;
      p3a::grid3 const cell_grid = dgt::generalize(block.cell_grid());
      double const detJ = block.cell_detJ();
      double const volume = dgt::get_volume(dim, block.domain().extents());
      View<double***> U = block.soln(0);
      auto f = [=] P3A_DEVICE (p3a::vector3<int> const& cell_ijk) {
        int const cell = cell_grid.index(cell_ijk);
        double val = 0.;
        for (int axis = 0; axis < dim; ++axis) {
          val += std::abs(U(cell, RH, 1 + axis)) * detJ;
        }
        return val;
      };
      double constexpr identity_value = p3a::zero_value<double>();
      auto constexpr binary_op = p3a::adder<double>();
      double const result = p3a::transform_reduce(
          p3a::execution::par, cell_grid, identity_value, binary_op, f);
      double const val = result / volume;
      if (val > tol) {
        marks[block.id()] = dgt::REFINE;
      }
    }
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
  } else if (amr == "sod") {
    refine_sod_initial(state.mesh);
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
  } else {
    throw std::runtime_error("invalid amr");
  }
}

}
