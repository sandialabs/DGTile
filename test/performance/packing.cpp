#include <chrono>

#include <dgt_cartesian.hpp>
#include <dgt_for_each.hpp>
#include <dgt_initialize.hpp>
#include <dgt_view.hpp>

#include <dgt_print.hpp> // debug

using namespace dgt;

static constexpr int num_blocks = 2;
static constexpr Grid3 cell_grid = {22, 22, 22};

struct Data
{
  int num_buffer_cells;
  View<real**> cell_values;
  View<real**> buffer_offsets;
  HostView<real**> buffer_offsets_h;
  View<Subgrid3*> owned_cells;
  HostView<Subgrid3*> owned_cells_h;
  View<real*> buffer;
};

static void setup_cell_values(Data& data)
{
  Kokkos::resize(data.cell_values, num_blocks, cell_grid.size());
  auto f = [=] DGT_DEVICE (int const b, Vec3<int> const& ijk) {
    int const cell = cell_grid.index(ijk);
    data.cell_values(b, cell) = cell;
  };
  for_each("setup_cell_values", num_blocks, cell_grid, f);
}

static void setup_buffer_offsets(Data& data)
{
  data.num_buffer_cells = 0;
  Kokkos::resize(data.buffer_offsets, num_blocks, offset_grid.size());
  Kokkos::resize(data.buffer_offsets_h, num_blocks, offset_grid.size());
  for (int b = 0; b < num_blocks; ++b) {
    auto f = [&] (Vec3<int> const& ijk) {
      if (ijk == Vec3<int>::zero()) return;
      int const idx = offset_grid.index(ijk);
      Vec3<std::int8_t> ijk_i8(ijk.x(), ijk.y(), ijk.z());
      Subgrid3 const owned_cells = get_cells(OWNED, cell_grid, ijk_i8);
      data.buffer_offsets_h(b, idx) = data.num_buffer_cells;
      data.num_buffer_cells += owned_cells.size();
    };
    seq_for_each(offset_grid, f);
  }
  Kokkos::deep_copy(data.buffer_offsets, data.buffer_offsets_h);
}

static void setup_owned_cells(Data& data)
{
  Kokkos::resize(data.owned_cells, offset_grid.size());
  Kokkos::resize(data.owned_cells_h, offset_grid.size());
  auto f = [&] (Vec3<int> const& ijk) {
    if (ijk == Vec3<int>::zero()) return;
    int const idx = offset_grid.index(ijk);
    Vec3<std::int8_t> ijk_i8(ijk.x(), ijk.y(), ijk.z());
    Subgrid3 const owned_cells = get_cells(OWNED, cell_grid, ijk_i8);
    data.owned_cells_h(idx) = owned_cells;
  };
  seq_for_each(offset_grid, f);
  Kokkos::deep_copy(data.owned_cells, data.owned_cells_h);
}

static void setup_buffer(Data& data)
{
  Kokkos::resize(data.buffer, data.num_buffer_cells);
}

static void setup_data(Data& data)
{
  setup_cell_values(data);
  setup_buffer_offsets(data);
  setup_owned_cells(data);
  setup_buffer(data);
}

static void pack_buffer_method_a(
    Data& data,
    int const block,
    Vec3<int> const& offset_ijk) {
  int const offset_idx = offset_grid.index(offset_ijk);
  Subgrid3 const owned_cells = data.owned_cells_h(offset_idx);
  int const buffer_offset = data.buffer_offsets_h(block, offset_idx);
  auto f = [=] DGT_DEVICE (Vec3<int> const& cell_ijk) {
    int const cell = cell_grid.index(cell_ijk);
    int const local_idx = owned_cells.index(cell_ijk);
    int const buffer_idx = buffer_offset + local_idx;
    data.buffer(buffer_idx) = data.cell_values(block, cell);
  };
  for_each("pack_method_a", owned_cells, f);
}

static std::int64_t pack_buffer_method_a(Data& data)
{
  auto const t0 = std::chrono::steady_clock::now();
  for (int b = 0; b < num_blocks; ++b) {
    auto f = [&] (Vec3<int> const& offset_ijk) {
      if (offset_ijk == Vec3<int>::zero()) return;
      pack_buffer_method_a(data, b, offset_ijk);
    };
    seq_for_each(offset_grid, f);
  }
  auto const t1 = std::chrono::steady_clock::now();
  auto const t = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();
  Kokkos::fence();
  printf(" > method a | %lld us\n", t);
  return t;
}

static std::int64_t pack_buffer_method_b(Data& data)
{
  std::vector<int> weights(offset_grid.size(), 0);
  for (int i = 0; i < offset_grid.size(); ++i) {
    weights[i] = data.owned_cells_h[i].size();
  }
  auto instances = Kokkos::Experimental::partition_space(
      Kokkos::DefaultExecutionSpace(),
      weights);
  return 1;
}

static std::int64_t pack_buffer_method_c()
{
  return 1;
}

static void do_packing_test()
{
  Data data;
  setup_data(data);
  printf("doing packing test\n");
  printf(" > num blocks: %d\n", num_blocks);
  printf(" > cell grid: [%d,%d,%d]\n",
      cell_grid.extents().x(),
      cell_grid.extents().y(),
      cell_grid.extents().z());
  std::int64_t a(0), b(0), c(0);
  a += pack_buffer_method_a(data);
  b += pack_buffer_method_b(data);
  c += pack_buffer_method_c();
  printf(" > total method a | %lld us\n", a);
  printf(" > total method b | %lld us\n", b);
  printf(" > total method c | %lld us\n", c);
}

int main(int argc, char** argv) {
  dgt::initialize(argc, argv);
  do_packing_test();
  dgt::finalize();
}
