#pragma once

#include <Kokkos_DualView.hpp>

#include "dgt_defines.hpp"
#include "dgt_dg.hpp"
#include "dgt_field.hpp"
#include "dgt_message.hpp"
#include "dgt_tree.hpp"
#include "dgt_subgrid3.hpp"

namespace dgt {

class Mesh;

struct Ghosting
{

  template <class T>
  using DualView = typename Kokkos::DualView<T>;

  using buffer_t = HostPinnedRightView<real***>;

  mpicpp::comm* m_comm = nullptr;
  Grid3 m_cell_grid = {0,0,0};

  int m_max_eqs = -1;
  int m_num_blocks = -1;
  int m_num_modes = -1;
  int m_num_msgs = -1;

  DualView<int*> m_block_offsets;
  DualView<int*> m_buffer_offsets;
  DualView<tree::Adjacent*> m_adjacencies;
  DualView<Subgrid3*> m_subgrids[DIRECTIONS];
  buffer_t m_buffers[DIRECTIONS];

  std::vector<Message<real>> m_messages[DIRECTIONS];

  void build(Mesh const& mesh);

  void begin_transfer(
      Field<real***> const& U,
      Basis<View> const& B,
      int const eq_start,
      int const eq_end);

  void end_transfer(
      Field<real***> const& U,
      Basis<View> const& B,
      int const eq_start,
      int const eq_end);


  void build_views(Mesh const& mesh);
  void build_messages(Mesh const& mesh);

  void set_message_size(int const neq);

  void pack(
      Field<real***> const& U,
      Basis<View> const& B,
      int const start_eq,
      int const end_eq);

  void post_messages();
  void wait_messages();

  void unpack(
      Field<real***> const& U,
      Basis<View> const& B,
      int const start_eq,
      int const end_eq);

};

}
