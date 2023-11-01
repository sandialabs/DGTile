#pragma once

#include <Kokkos_DualView.hpp>

#include "dgt_defines.hpp"
#include "dgt_message.hpp"
#include "dgt_subgrid3.hpp"
#include "dgt_tree.hpp"
#include "dgt_view.hpp"

namespace dgt {

class Mesh;

class Ghosting
{

  private:

    template <class T>
    using DualView = typename Kokkos::DualView<T>;

    using buffer_t = HostPinnedView<real***>;

    int m_max_eqs = -1;
    int m_num_modes = -1;

    DualView<int*> m_msgs_per_block;
    DualView<int*> m_buffer_offsets;
    DualView<tree::Adjacent*> m_adjacencies;
    DualView<Subgrid3*> m_subgrids[DIRECTIONS];
    buffer_t m_buffers[DIRECTIONS];

    std::vector<Message<real>> m_messages[DIRECTIONS];

  public:

    void build(Mesh const& mesh);

  private:

    void build_views(Mesh const& mesh);
    void build_messages(Mesh const& mesh);

};

}
