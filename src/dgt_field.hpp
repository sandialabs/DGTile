#pragma once

#include "dgt_defines.hpp"
#include "dgt_view.hpp"

namespace dgt {

class ModalField
{

  public:

    using view_t = View<real***>;
    using uview_t = UnmanagedView<real***>;
    using storage_t = std::vector<view_t>;
    using accessor_t = HostPinnedView<uview_t*>;

  private:

    storage_t m_storage;
    accessor_t m_accessor;

  public:

    ModalField(
        std::string const& name,
        int const num_blocks,
        int const num_cells,
        int const num_eqs,
        int const num_modes);

    ~ModalField();

    accessor_t get() { return m_accessor; }

};

}
