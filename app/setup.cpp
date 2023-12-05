#include "app.hpp"

namespace app {

void setup(mpicpp::comm* comm, Input const& in, State& state)
{
  int const p = in.basis.polynomial_order;
  int const q = in.basis.quadrature_rule;
  bool const tensor = in.basis.tensor_product;
  state.mesh.set_comm(comm);
  state.mesh.set_cell_grid(in.mesh.cell_grid);
  state.mesh.set_periodic(in.mesh.periodic);
  state.mesh.set_basis(p, q, tensor);
  state.integrator = create_integrator(in.time.integrator);
  //TODO: set up some integrators/physics here
}

}
