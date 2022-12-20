#pragma once

#include "p3a_grid3.hpp"

namespace dgt {

using namespace p3a;

struct Point;

[[nodiscard]] int get_dim(vector3<int> const& v);
[[nodiscard]] int get_dim(grid3 const& g);

[[nodiscard]] grid3 generalize(grid3 const& g);
[[nodiscard]] subgrid3 generalize(subgrid3 const& s);

[[nodiscard]] grid3 get_child_grid(int dim);

[[nodiscard]] grid3 get_side_grid(grid3 const& cells, int axis);
[[nodiscard]] subgrid3 get_intr_sides(grid3 const& cells, int axis);
[[nodiscard]] subgrid3 get_adj_sides(grid3 const& cells, int axis, int dir);
[[nodiscard]] subgrid3 get_adj_cells(grid3 const& cells, int axis, int dir);

[[nodiscard]] grid3 get_viz_cell_grid(grid3 const& cells, int p);
[[nodiscard]] grid3 get_viz_point_grid(grid3 const& cells, int p);

[[nodiscard]] bool contains(grid3 const& g, subgrid3 const& s);
[[nodiscard]] bool contains(int depth, subgrid3 const& s, Point const& pt);

[[nodiscard]] int get_num_local(int ntotal, int nparts, int part);
[[nodiscard]] int get_local_offset(int ntotal, int nparts, int part);

}
