#include <dgt_cartesian.hpp>

#include "example.hpp"

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

void setup(State& state, mpicpp::comm* comm, Input const& in)
{
  int const dim = infer_dimension(in.mesh.cell_grid);
  Basis<View> basis = build_basis<View>(
      dim,
      in.basis.polynomial_order,
      in.basis.quadrature_rule,
      in.basis.tensor_product);
  state.eqs = Equations(in.num_materials);
  state.mesh.set_comm(comm);
  state.mesh.set_domain(in.mesh.domain);
  state.mesh.set_cell_grid(in.mesh.cell_grid);
  state.mesh.set_periodic(in.mesh.periodic);
  state.mesh.set_basis(basis);
  state.mesh.init(in.mesh.block_grid);
  int const nstored = get_num_stored_solutions(basis.p);
  int const neqs = state.eqs.num_eqs();
  state.mesh.add_modal("hydro", nstored, neqs);
}

}
