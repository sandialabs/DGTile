#pragma once

#include <string>

#include "dgt_defines.hpp"
#include "dgt_grid3.hpp"
#include "dgt_view.hpp"

namespace dgt {
namespace field {

template <class T>
class Base
{
  public:
    using view_t = View<T>;
    using uview_t = UnmanagedView<T>;
    using storage_t = std::vector<view_t>;
    using accessor_t = HostPinnedView<uview_t*>;
  protected:
    std::string m_name;
    storage_t m_storage;
    accessor_t m_accessor;
  public:
    std::string name() const { return m_name; }
    accessor_t get() { return m_accessor; }
    view_t get_view(int const block) { return m_storage[block]; }
};

template <class T>
class AxisBase
{
  public:
    using view_t = View<T>;
    using uview_t = UnmanagedView<T>;
    using storage_t = std::vector<view_t>;
    using accessor_t = HostPinnedView<uview_t*>;
  protected:
    std::string m_name;
    Vec3<storage_t> m_storage;
    Vec3<accessor_t> m_accessor;
  public:
    std::string name() const { return m_name; }
    Vec3<accessor_t> get() { return m_accessor; }
    view_t get_view(int const axis, int const block) { return m_storage[axis][block]; }
};

class Modal : public Base<real***>
{
  public:
    Modal(
        std::string const& name,
        Grid3 const& cell_grid,
        int const num_blocks,
        int const num_comps,
        int const num_modes);
};

class Cell : public Base<real**>
{
  public:
    Cell(
        std::string const& name,
        Grid3 const& cell_grid,
        int const num_blocks,
        int const num_comps);
};

class CellPoints : public Base<real***>
{
  public:
    CellPoints(
        std::string const& name,
        Grid3 const& cell_grid,
        int const num_blocks,
        int const num_points,
        int const num_comps);
};

class Face : public AxisBase<real**>
{
  public:
    Face(
        std::string const& name,
        Grid3 const& cell_grid,
        int const num_blocks,
        int const num_comps);
};

class FacePoints : public AxisBase<real***>
{
  public:
    FacePoints(
        std::string const& name,
        Grid3 const& cell_grid,
        int const num_blocks,
        int const num_points,
        int const num_comps);
};

}
}
