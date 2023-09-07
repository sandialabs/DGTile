#pragma once

#include "dgt_defines.hpp"
#include "dgt_dg.hpp"
#include "dgt_view.hpp"

namespace dgt {

class DGField 
{

  public:

    using view_t = View<real***>;
    using uview_t = UnmanagedView<real***>;
    using storage_t = std::vector<view_t>;
    using accessor_t = HostPinnedView<uview_t*>;

  private:

    std::string m_name;
    BasisInfo m_basis;

    storage_t m_storage;
    accessor_t m_accessor;

  public:

    DGField(
        std::string const& name,
        int const num_blocks,
        int const num_cells,
        int const num_eqs,
        BasisInfo const& basis);

    ~DGField();

    std::string name() const { return m_name; }
    BasisInfo basis() const { return m_basis; }

    accessor_t get() { return m_accessor; }
    view_t get(int const idx) { return m_storage[idx]; }

};

}
