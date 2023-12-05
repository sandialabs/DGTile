#include "app.hpp"
#include "hydro.hpp"

namespace app {

static void describe_integrator(State const& state)
{
  if (state.mesh.comm()->rank()) return;
  printf("integrator:\n");
  printf("> %s\n", state.integrator->name().c_str());
  printf("> stages: %d\n", state.integrator->num_stages());
  printf("> required storage: %d\n", state.integrator->required_containers());
}

void setup(mpicpp::comm* comm, Input const& in, State& state)
{
  int const p = in.basis.polynomial_order;
  int const q = in.basis.quadrature_rule;
  bool const tensor = in.basis.tensor_product;
  state.mesh.set_comm(comm);
  state.mesh.set_domain(in.mesh.domain);
  state.mesh.set_cell_grid(in.mesh.cell_grid);
  state.mesh.set_periodic(in.mesh.periodic);
  state.mesh.set_basis(p, q, tensor);
  state.integrator = create_integrator(in.time.integrator);
  state.physics.push_back(std::make_unique<Hydro>(&state, &in));
  state.mesh.initialize(in.mesh.block_grid);
  describe_integrator(state);
  state.mesh.print_stats();
}

}
