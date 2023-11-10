#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint.hpp>
#include <conduit/conduit_relay.hpp>
#include <conduit/conduit_relay_io_silo.hpp>

#include <fmt/core.h>

#include "dgt_cartesian.hpp"
#include "dgt_conduit.hpp"
#include "dgt_mesh.hpp"

#include "dgt_print.hpp"

namespace dgt {
namespace silo {

void dummy(Mesh const& mesh)
{
  conduit::Node node;
  int const dim = mesh.dim();
  int const num_blocks = mesh.num_owned_blocks();
  Grid3 const cell_grid = mesh.cell_grid();
  Subgrid3 const owned_cells = get_owned_cells(cell_grid);
  Vec3<int> const ncells = owned_cells.extents();
  std::string const dim_names[DIMENSIONS] = {"i", "j", "k"};
  std::string const origin_names[DIMENSIONS] = {"x", "y", "z"};
  std::string const spacing_names[DIMENSIONS] = {"dx", "dy", "dz"};
  for (int block = 0; block < num_blocks; ++block) {
    Box3<real> const domain = mesh.block_info_h().domains[block];
    Vec3<real> const dx = mesh.block_info_h().cell_dxs[block];
    Vec3<real> const o = domain.lower() + dx;
    std::string const block_name = fmt::format("block_{}", block);
    conduit::Node& b = node[block_name];
//    b["state/cycle"] = 0;
    b["state/domain_id"] = block;
    b["coordsets/coords/type"] = "uniform";
    b["topologies/topo/type"] = "uniform";
    b["topologies/topo/coordset"] = "coords";
    for (int axis = 0; axis < dim; ++axis) {
      b["coordsets/coords/dims/" + dim_names[axis]] = ncells[axis] + 1;
      b["coordsets/coords/origin/" + origin_names[axis]] = o[axis];
      b["coordsets/coords/spacing/" + spacing_names[axis]] = dx[axis];
    }
  }
  conduit::Node debug;
  if (!conduit::blueprint::mesh::verify(node, debug)) {
    std::string const msg = fmt::format(
        "dgt::silo::dummy(): {}\n", debug.to_yaml());
    throw std::runtime_error(msg);
  }
  conduit::Node options;
  options["number_of_files"] = 1;
  conduit::relay::io::silo::save_mesh(node, "debug", options);
}

}
}
