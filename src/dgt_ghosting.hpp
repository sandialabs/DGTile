#pragma once

#include <Kokkos_DualView.hpp>

#include <mpicpp.hpp>

#include "dgt_field.hpp"
#include "dgt_subgrid3.hpp"
#include "dgt_tree.hpp"
#include "dgt_view.hpp"

namespace dgt {

class Mesh;

struct Message
{
  int tag = -1;
  int rank = -1;
  int size = -1;
  real* data = nullptr;
  mpicpp::request req;
  void send(mpicpp::comm* c) { req = c->isend(data, size, rank, tag); }
  void recv(mpicpp::comm* c) { req = c->irecv(data, size, rank, tag); }
  void wait() { req.wait(); }
};

class Ghosting
{

  private:

    template <class T>
    using DualView = typename Kokkos::DualView<T>;

    Grid3 m_cell_grid;

    int m_num_messages = -1;
    int m_num_cells = -1;

    DualView<int*> m_msgs_per_block;
    DualView<Subgrid3*> m_subgrids[DIRECTIONS];
    DualView<int*> m_buffer_offsets[DIRECTIONS];
    std::vector<Message> m_messages[DIRECTIONS];
    HostPinnedRightView<real***> m_buffer[DIRECTIONS];

  public:

    void build(dgt::Mesh const& mesh);

  private:

    void build_helper_views(dgt::Mesh const& mesh);
    void build_buffers(dgt::Mesh const& mesh);
    void build_messages(dgt::Mesh const& mesh);


};

}
