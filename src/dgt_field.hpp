#pragma once

#include "dgt_defines.hpp"
#include "dgt_view.hpp"

namespace dgt {

class ModalField
{

  private:

    using view_t = View<real***>;
    using accessor_t = HostPinnedView<view_t*>;
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
