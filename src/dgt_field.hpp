#pragma once

#include "dgt_defines.hpp"
#include "dgt_view.hpp"

namespace dgt {

class Field
{

  public:

    using view_t = View<real***>;
    using uview_t = UnmanagedView<real***>;
    using storage_t = std::vector<view_t>;
    using accessor_t = HostPinnedView<uview_t*>;

  public:

    enum {MODAL, FLUX, RESIDUAL, CELLS, FACES, CELL_POINTS, FACE_POINTS, KINDS};

  private:

    std::string m_name;

    storage_t m_storage;
    accessor_t m_accessor;

  public:

    Field(
        std::string const& name,
        int const num_blocks,
        int const extent0,
        int const extent1,
        int const extent2);

    std::string name() const { return m_name; }

    accessor_t get() { return m_accessor; }
    view_t get(int const idx) { return m_storage[idx]; }

};

}
