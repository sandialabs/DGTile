#include <fmt/core.h>

#include "dgt_cartesian.hpp"
#include "dgt_field.hpp"

namespace dgt {
namespace field {

template <class BaseT>
void allocate(
    typename BaseT::storage_t& storage,
    typename BaseT::accessor_t& accessor,
    std::string const& name,
    int const view_dim,
    int const num_blocks,
    int const n0,
    int const n1,
    int const n2)
{
  using Kokkos::view_alloc;
  using Kokkos::WithoutInitializing;
  storage.resize(num_blocks);
  accessor = typename BaseT::accessor_t(name, num_blocks);
  for (int block = 0; block < num_blocks; ++block) {
    auto const block_name = fmt::format("{}[{}]", name, block);
    auto const alloc = view_alloc(block_name, WithoutInitializing);
    if (view_dim == 2) {
      storage[block] = typename BaseT::view_t(alloc, n0, n1);
      accessor[block] = typename BaseT::uview_t(storage[block].data(), n0, n1);
    }
    if (view_dim == 3) {
      storage[block] = typename BaseT::view_t(alloc, n0, n1, n2);
      accessor[block] = typename BaseT::uview_t(storage[block].data(), n0, n1, n2);
    }
  }
}

Modal::Modal(
    std::string const& name,
    Grid3 const& cell_grid,
    int const num_blocks,
    int const num_comps,
    int const num_modes)
{
  int const view_dim = 3;
  int const num_cells = get_num_cells(cell_grid);
  allocate<field::Modal>(m_storage, m_accessor,
      name, view_dim, num_blocks, num_cells, num_comps, num_modes);
}

Cell::Cell(
    std::string const& name,
    Grid3 const& cell_grid,
    int const num_blocks,
    int const num_comps)
{
  int const view_dim = 2;
  int const num_cells = get_num_cells(cell_grid);
  allocate<field::Cell>(m_storage, m_accessor,
      name, view_dim, num_blocks, num_cells, num_comps, 0);
}

CellPoints::CellPoints(
    std::string const& name,
    Grid3 const& cell_grid,
    int const num_blocks,
    int const num_points,
    int const num_comps)
{
  int const view_dim = 3;
  int const num_cells = get_num_cells(cell_grid);
  allocate<field::CellPoints>(m_storage, m_accessor,
      name, view_dim, num_blocks, num_cells, num_points, num_comps);
}

Face::Face(
    std::string const& name,
    Grid3 const& cell_grid,
    int const num_blocks,
    int const num_comps)
{
  int const view_dim = 2;
  int const dim = infer_dimension(cell_grid);
  for (int axis = 0; axis < dim; ++axis) {
    auto const axis_name = fmt::format("{}{}", name, get_axis_name(axis));
    int const num_faces = get_num_faces(cell_grid, axis);
    allocate<field::Face>(m_storage[axis], m_accessor[axis],
        axis_name, view_dim, num_blocks, num_faces, num_comps, 0);
  }
}

FacePoints::FacePoints(
    std::string const& name,
    Grid3 const& cell_grid,
    int const num_blocks,
    int const num_points,
    int const num_comps)
{
  int const view_dim = 3;
  int const dim = infer_dimension(cell_grid);
  for (int axis = 0; axis < dim; ++axis) {
    auto const axis_name = fmt::format("{}{}", name, get_axis_name(axis));
    int const num_faces = get_num_faces(cell_grid, axis);
    allocate<field::FacePoints>(m_storage[axis], m_accessor[axis],
        axis_name, view_dim, num_blocks, num_faces, num_points, num_comps);
  }
}

}
}
