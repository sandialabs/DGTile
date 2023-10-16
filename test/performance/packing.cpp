#include <chrono>

#include <dgt_cartesian.hpp>
#include <dgt_for_each.hpp>
#include <dgt_initialize.hpp>
#include <dgt_view.hpp>

using namespace dgt;
using namespace std::chrono;

struct Data
{
  int num_blocks;
  Grid3 cell_grid;
  int num_messages;
  int num_message_cells;
  View<real**> input_values;
  View<int*> message_blocks;
  HostView<int*> message_blocks_h;
  View<int*> message_offsets;
  HostView<int*> message_offsets_h;
  View<Subgrid3*> message_subgrids;
  HostView<Subgrid3*> message_subgrids_h;
  View<real*> packed_values;
};

static void setup_input_values(Data& data)
{
  Kokkos::resize(data.input_values, data.num_blocks, data.cell_grid.size());
  Grid3 const cgrid = data.cell_grid;
  auto input_values = data.input_values;
  auto f = [=] DGT_DEVICE (int const b, Vec3<int> const& ijk) DGT_ALWAYS_INLINE {
    int const cell = cgrid.index(ijk);
    input_values(b, cell) = cell;
  };
  for_each("setup_cell_values", data.num_blocks, data.cell_grid, f);
  Kokkos::fence();
}

static void setup_messages(Data& data)
{
  // this would be some separate counting function with AMR
  data.num_messages = data.num_blocks * (meta_grid.size()-1);
  data.num_message_cells = 0;
  Kokkos::resize(data.message_blocks, data.num_messages);
  Kokkos::resize(data.message_blocks_h, data.num_messages);
  Kokkos::resize(data.message_offsets, data.num_messages);
  Kokkos::resize(data.message_offsets_h, data.num_messages);
  Kokkos::resize(data.message_subgrids, data.num_messages);
  Kokkos::resize(data.message_subgrids_h, data.num_messages);
  int msg_idx = 0;
  for (int b = 0; b < data.num_blocks; ++b) {
    auto f = [&] (Vec3<int> const& ijk) {
      if (ijk == Vec3<int>::zero()) return;
      Vec3<std::int8_t> ijk_i8(ijk.x(), ijk.y(), ijk.z());
      Subgrid3 const owned_cells = get_cells(OWNED, data.cell_grid, ijk_i8);
      data.message_blocks_h[msg_idx] = b;
      data.message_subgrids_h[msg_idx] = owned_cells;
      data.message_offsets_h[msg_idx] = data.num_message_cells;
      data.num_message_cells += owned_cells.size();
      msg_idx++;
    };
    seq_for_each(meta_grid, f);
  }
  Kokkos::deep_copy(data.message_blocks, data.message_blocks_h);
  Kokkos::deep_copy(data.message_offsets, data.message_offsets_h);
  Kokkos::deep_copy(data.message_subgrids, data.message_subgrids_h);
  Kokkos::fence();
}

static void setup_packed_values(Data& data)
{
  Kokkos::resize(data.packed_values, data.num_message_cells);
}

static void setup_data(Data& data)
{
  setup_input_values(data);
  setup_messages(data);
  setup_packed_values(data);
}

// launch an individual kernel for each block's border
static std::int64_t pack_values_method_a(Data& data)
{
  auto const t0 = steady_clock::now();
  Grid3 const cell_grid = data.cell_grid;
  auto input_values = data.input_values;
  auto packed_values = data.packed_values;
  int msg_idx = 0;
  for (int b = 0; b < data.num_blocks; ++b) {
    auto f1 = [&] (Vec3<int> const& offset_ijk) DGT_ALWAYS_INLINE {
      if (offset_ijk == Vec3<int>::zero()) return;
      int const packed_offset = data.message_offsets_h[msg_idx];
      Subgrid3 const owned_cells = data.message_subgrids_h[msg_idx];
      auto f2 = [=] DGT_DEVICE (Vec3<int> const& cell_ijk) DGT_ALWAYS_INLINE {
        int const cell = cell_grid.index(cell_ijk);
        int const local_idx = owned_cells.index(cell_ijk);
        int const packed_idx = packed_offset + local_idx;
        packed_values[packed_idx] = input_values(b, cell);
      };
      for_each("", owned_cells, f2);
      msg_idx++;
    };
    seq_for_each(meta_grid, f1);
  }
  Kokkos::fence();
  auto const t1 = steady_clock::now();
  auto const t = duration_cast<microseconds>(t1-t0).count();
  return t;
}

DGT_METHOD inline bool contains(
    Subgrid3 const& s,
    Vec3<int> const& ijk)
{
  return
    ijk.x() >= s.lower().x() &&
    ijk.y() >= s.lower().y() &&
    ijk.z() >= s.lower().z() &&
    ijk.x() < s.upper().x() &&
    ijk.y() < s.upper().y() &&
    ijk.z() < s.upper().z();
}

// loop over all cells in the mesh and check if a cell is in a border
static std::int64_t pack_values_method_b(Data& data)
{
  auto const t0 = steady_clock::now();
  int const nmsg_per_block = meta_grid.size()-1;
  Grid3 const cell_grid = data.cell_grid;
  auto message_offsets = data.message_offsets;
  auto message_subgrids = data.message_subgrids;
  auto input_values = data.input_values;
  auto packed_values = data.packed_values;
  auto f = [=] DGT_DEVICE (int const block, Vec3<int> const& cell_ijk) DGT_ALWAYS_INLINE
  {
    int const cell = cell_grid.index(cell_ijk);
    for (int i = 0; i < nmsg_per_block; ++i) { // this would be more complicated w/ AMR
      int const msg_idx = i + block*nmsg_per_block;
      Subgrid3 const subgrid = message_subgrids[msg_idx];
      if (!contains(subgrid, cell_ijk)) continue;
      int const packed_offset = message_offsets[msg_idx];
      int const local_idx = subgrid.index(cell_ijk);
      int const packed_idx = packed_offset + local_idx;
      packed_values[packed_idx] = input_values(block, cell);
    }
  };
  for_each("", data.num_blocks, data.cell_grid, f);
  Kokkos::fence();
  auto const t1 = steady_clock::now();
  auto const t = duration_cast<microseconds>(t1-t0).count();
  return t;
}

// use fancy kokkos teams
static std::int64_t pack_values_method_c(Data& data)
{
  auto const t0 = steady_clock::now();
  auto policy = Kokkos::TeamPolicy<>(data.num_messages, Kokkos::AUTO());
  using member_type = typename decltype(policy)::member_type;
  auto cell_grid = data.cell_grid;
  auto message_blocks = data.message_blocks;
  auto message_offsets = data.message_offsets;
  auto message_subgrids = data.message_subgrids;
  auto input_values = data.input_values;
  auto packed_values = data.packed_values;
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (member_type member) {
    int msg_idx = member.league_rank();
    int block_idx = message_blocks[msg_idx];
    int packed_offset = message_offsets[msg_idx];
    Subgrid3 const subgrid = message_subgrids[msg_idx];
    Vec3<int> const ncells = subgrid.extents();
    Grid3 const subgrid_indexer(subgrid.extents());
    auto range = Kokkos::TeamThreadMDRange<Kokkos::Rank<3>, member_type>(member, ncells.x(), ncells.y(), ncells.z());
    Kokkos::parallel_for(range, [=](int i, int j, int k) {
      Vec3<int> cell_ijk = subgrid.lower() + Vec3<int>(i,j,k);
      int const packed_idx = packed_offset + subgrid.index(cell_ijk);
      int const cell_idx = cell_grid.index(cell_ijk);
      packed_values[packed_idx] = input_values(block_idx, cell_idx);
    });
  });
  Kokkos::fence();
  auto const t1 = steady_clock::now();
  auto const t = duration_cast<microseconds>(t1-t0).count();
  return t;
}

static void do_packing_test(Data& data)
{
  setup_data(data);
  printf("doing packing test\n");
  printf(" > num blocks: %d\n", data.num_blocks);
  printf(" > cell grid: [%d,%d,%d]\n", data.cell_grid.extents().x(), data.cell_grid.extents().y(), data.cell_grid.extents().z());
  printf(" > num message cells: %d\n", data.num_message_cells);
  std::int64_t a(0), b(0), c(0);
  for (int i = 1; i <= 20; ++i) {
    Kokkos::deep_copy(data.packed_values, -1.);
    Kokkos::fence();
    a += pack_values_method_a(data);
    Kokkos::deep_copy(data.packed_values, -1.);
    Kokkos::fence();
    b += pack_values_method_b(data);
    Kokkos::deep_copy(data.packed_values, -1.);
    Kokkos::fence();
    c += pack_values_method_c(data);
  }
  printf(" > total method a | %lld us\n", a);
  printf(" > total mehtod b | %lld us\n", b);
  printf(" > total method c | %lld us\n", c);
}

int main(int argc, char** argv) {
  dgt::initialize(argc, argv);
  int const nblocks = std::stoi(argv[1]);
  int const ncells = std::stoi(argv[2]);
  {
    Data data;
    data.num_blocks = nblocks;
    data.cell_grid = Grid3(ncells, ncells, ncells);
    do_packing_test(data);
  }
  dgt::finalize();
}
