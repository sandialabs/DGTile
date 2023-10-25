#include <dgt_cartesian.hpp>
#include <dgt_dg.hpp>
#include <dgt_for_each.hpp>

#include "example.hpp"

namespace example {

static int get_num_stored_solutions(int const p)
{
  static constexpr int table[max_polynomial_order+1] = {1, 2, 2, 3};
  return table[p];
}

static void apply_initial_conditions(State& state, Input const& in)
{
  static constexpr int CELL = basis_locations::CELL;
  Mesh& mesh = state.mesh;
  int const nblocks = mesh.num_owned_blocks();
  Field<real***> U_field = mesh.get_solution("hydro", 0);
  Grid3 const cell_grid = mesh.cell_grid();
  int const ncells = generalize(mesh.dim(), cell_grid).size();
  EoS const& eos = state.eos;
  inputs::InitialConditions const& ics = in.ics;
  Basis<HostView> const& B = mesh.basis_h();
  BlockInfo<HostView> const& info = mesh.block_info_h();
  HostView<real***> U_host("U_ic", ncells, NEQ, B.num_modes);
  for (int block = 0; block < nblocks; ++block) {
    Kokkos::deep_copy(U_host, 0.);
    Vec3<real> const origin = info.domains[block].lower();
    Vec3<real> const dx = info.cell_dxs[block];
    auto functor = [&] (Vec3<int> const& cell_ijk) {
      int const cell = cell_grid.index(cell_ijk);
      for (int pt = 0; pt < B.num_cell_pts; ++pt) {
        real const wt = B.cell_weights(pt);
        Vec3<real> const xi = get_point(B, CELL, pt);
        Vec3<real> const x = map_to_physical(cell_ijk, origin, dx, xi);
        real const rho = ics.density->operator()(x);
        real const p = ics.pressure->operator()(x);
        Vec3<real> const v = ics.velocity->operator()(x);
        Vec3<real> const mmtm = rho * v;
        real const e = eos.e_from_rho_p(rho, p);
        real const half_v2 = 0.5 * dot(v, v);
        real const En = rho*e + rho*half_v2;
        for (int mode = 0; mode < B.num_modes; ++mode) {
          real const phi = B.modes[CELL].phis(pt, mode);
          real const M_inv = 1. / B.mass(mode);
          U_host(cell, DENS, mode) += rho * phi * wt * M_inv;
          U_host(cell, MMTM + X, mode) += mmtm.x() * phi * wt * M_inv;
          U_host(cell, MMTM + Y, mode) += mmtm.y() * phi * wt * M_inv;
          U_host(cell, MMTM + Z, mode) += mmtm.z() * phi * wt * M_inv;
          U_host(cell, ENER, mode) += En * phi * wt * M_inv;
        }
      }
    };
    seq_for_each(cell_grid, functor);
    Kokkos::deep_copy(U_field.get_view(block), U_host);
  }
}

void setup(State& state, mpicpp::comm* comm, Input const& in)
{
  int const p = in.basis.polynomial_order;
  int const q = in.basis.quadrature_rule;
  bool const tensor = in.basis.tensor_product;
  state.mesh.set_comm(comm);
  state.mesh.set_domain(in.mesh.domain);
  state.mesh.set_cell_grid(in.mesh.cell_grid);
  state.mesh.set_periodic(in.mesh.periodic);
  state.mesh.set_basis(p, q, tensor);
  int const nstored = get_num_stored_solutions(p);
  state.mesh.add_modal({"hydro", nstored, NEQ, true});
  state.mesh.initialize(in.mesh.block_grid);
  state.mesh.print_stats();
  state.eos = EoS(in.gamma);
  apply_initial_conditions(state, in);

  // debug
  write_out(in, state, 0);

}

}
