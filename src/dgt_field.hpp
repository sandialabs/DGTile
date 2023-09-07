#pragma once

#include <string>

#include "dgt_defines.hpp"
#include "dgt_grid3.hpp"
#include "dgt_view.hpp"

namespace dgt {
namespace field {

template <class T>
class CellBase
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

  protected:

    void allocate(
        int const num_blocks,
        int const extent0,
        int const extent1,
        int const extent2);

};

class Modal : public CellBase<real***>
{
  public:
    Modal(
        std::string const& name,
        Grid3 const& cell_grid,
        int const num_blocks,
        int const num_comps,
        int const num_modes);
};

class Cell : public CellBase<real**>
{
  public:
    Cell(
        std::string const& name,
        Grid3 const& cell_grid,
        int const num_blocks,
        int const num_comps);
};

class CellPoints : public CellBase<real***>
{
  public:
    CellPoints(
        std::string const& name,
        Grid3 const& cell_grid,
        int const num_blocks,
        int const num_comps);
};

}
}
