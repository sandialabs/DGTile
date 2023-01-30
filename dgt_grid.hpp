#pragma once

#include "p3a_grid3.hpp"

namespace dgt {

struct Point;

[[nodiscard]] int get_dim(p3a::vector3<int> const& v);
[[nodiscard]] int get_dim(p3a::grid3 const& g);

[[nodiscard]] p3a::grid3 generalize(p3a::grid3 const& g);
[[nodiscard]] p3a::subgrid3 generalize(p3a::subgrid3 const& s);

[[nodiscard]] p3a::grid3 get_child_grid(int dim);

[[nodiscard]] p3a::grid3 get_node_grid(p3a::grid3 const& cells);
[[nodiscard]] p3a::grid3 get_side_grid(p3a::grid3 const& cells, int axis);
[[nodiscard]] p3a::subgrid3 get_intr_sides(p3a::grid3 const& cells, int axis);
[[nodiscard]] p3a::subgrid3 get_adj_sides(p3a::grid3 const& cells, int axis, int dir);
[[nodiscard]] p3a::subgrid3 get_adj_cells(p3a::grid3 const& cells, int axis, int dir);

[[nodiscard]] bool contains(p3a::grid3 const& g, p3a::subgrid3 const& s);
[[nodiscard]] bool contains(int depth, p3a::subgrid3 const& s, Point const& pt);

[[nodiscard]] int get_num_local(int ntotal, int nparts, int part);
[[nodiscard]] int get_local_offset(int ntotal, int nparts, int part);

}
