#include <dgt_cartesian.hpp>
#include <dgt_for_each.hpp>

#include "example.hpp"

#include <dgt_print.hpp> // debug

namespace example {

Equations::Equations(int const num_mats)
{
  std::fill_n(offsets, NVAR, -1);
  offsets[RHO] = 0;
  offsets[MMTM] = num_mats;
  offsets[ENER] = num_mats + DIMENSIONS;
}

static int get_num_stored_solutions(int const p)
{
  static constexpr int table[max_polynomial_order+1] = {1, 2, 2, 3};
  return table[p];
}

static void apply_initial_conditions(State& state, Input const& in)
{
  Mesh& mesh = state.mesh;
  int const nmats = in.num_materials;
  int const nblocks = mesh.num_owned_blocks();
  Equations const& eqs = state.eqs;
  Field<real***> U = mesh.get_solution("hydro", 0);
  Grid3 const cell_grid = mesh.cell_grid();
  int const ncells = generalize(mesh.dim(), cell_grid).size();
  Basis<HostView> const& B = mesh.basis_h();
  HostView<real***> U_host("U_ic", ncells, eqs.num_eqs(), B.num_modes);
  inputs::InitialConditions const& ics = in.ics;
  inputs::Materials const mats = in.materials;
  for (int block = 0; block < nblocks ; ++block) {
    Kokkos::deep_copy(U_host, 0.);
    auto functor = [&] (Vec3<int> const& cell_ijk) {
      int const cell = cell_grid.index(cell_ijk);
      for (int pt = 0; pt < B.num_cell_pts; ++pt) {
        Vec3<real> const xi = get_point(B, basis_locations::CELL, pt);
        Vec3<real> v = ics.velocity->operator()(xi);
        real const half_v2 = 0.5*dot(v, v);
        real rho_bulk = 0.;
        for (int mat = 0; mat < nmats; ++mat) {
          real const rho = ics.densities[mat]->operator()(xi);
        }
      }
    };
    seq_for_each(cell_grid, functor);
    Kokkos::deep_copy(U.get_view(block), U_host);
  }
}

void setup(State& state, mpicpp::comm* comm, Input const& in)
{
  int const dim = infer_dimension(in.mesh.cell_grid);
  int const p = in.basis.polynomial_order;
  int const q = in.basis.quadrature_rule;
  bool const tensor = in.basis.tensor_product;
  state.eqs = Equations(in.num_materials);
  state.mesh.set_comm(comm);
  state.mesh.set_domain(in.mesh.domain);
  state.mesh.set_cell_grid(in.mesh.cell_grid);
  state.mesh.set_periodic(in.mesh.periodic);
  state.mesh.set_basis(p, q, tensor);
  state.mesh.init(in.mesh.block_grid);
  state.mesh.print_stats();
  int const nstored = get_num_stored_solutions(p);
  int const neqs = state.eqs.num_eqs();
  state.mesh.add_modal("hydro", nstored, neqs);
  apply_initial_conditions(state, in);
}

}
