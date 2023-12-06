#include <filesystem>

#include <dgt_cartesian.hpp>
#include <dgt_for_each.hpp>
#include <dgt_initialize.hpp>
#include <dgt_mesh.hpp>
#include <dgt_stl.hpp>
#include <dgt_vtk.hpp>

#include <dgt_print.hpp> // debug

using namespace dgt;
using namespace dgt::stl;

static void print_usage(std::string const& exe)
{
  printf("usage: %s <geom.stl> <cells_per_axis>\n", exe.c_str());
}

static Mesh make_mesh(
    mpicpp::comm* comm,
    Box3<real> const& domain,
    int const cells)
{
  Mesh mesh;
  mesh.set_comm(comm);
  mesh.set_domain(domain);
  mesh.set_cell_grid({cells, cells, cells});
  mesh.set_periodic({false, false, false});
  mesh.set_basis(1, 0, true);
  mesh.add_modal({"vfs", 1, 1, true});
  mesh.initialize(Grid3(1,1,1));
  return mesh;
}

dgt::vtk::VtkView<real> get_vfs(
    Mesh const& mesh,
    Field<real**> const& vf_field,
    int const block)
{
  Basis<View> const& B = mesh.basis();
  Grid3 const cell_grid = mesh.cell_grid();
  Grid3 const inner_grid = tensor_bounds(B.dim, B.q-1);
  Grid3 const ginner_grid = generalize(B.dim, inner_grid);
  Grid3 const viz_cell_grid = vtk::get_viz_cell_grid(cell_grid, B.q);
  Grid3 const gviz_cell_grid = generalize(B.dim, viz_cell_grid);
  Subgrid3 const owned_cells = get_owned_cells(cell_grid);
  Vec3<int> const ghost_offset = owned_cells.lower();
  auto const vfs = vf_field.get();
  dgt::vtk::VtkView<real> vf_viz;
  Kokkos::resize(vf_viz, gviz_cell_grid.size(), 1);
  auto functor = [=] DGT_HOST_DEVICE (Vec3<int> const& cell_ijk) {
    int const cell = cell_grid.index(cell_ijk);
    Vec3<int> const owned_ijk = cell_ijk - ghost_offset;
    inner_for_each(ginner_grid,
    [&] (Vec3<int> const& inner_ijk) DGT_ALWAYS_INLINE {
      int const pt = ginner_grid.index(inner_ijk);
      Vec3<int> const viz_cell_ijk = (B.q * owned_ijk) + inner_ijk;
      int const viz_cell = gviz_cell_grid.index(viz_cell_ijk);
      real const val = vfs[block](cell, pt);
      vf_viz.d_view(viz_cell, 0) = val;
    });
  };
  for_each("vtk::get_variable", owned_cells, functor);
  return vf_viz;
}

static void write_mesh(
    std::filesystem::path const& path,
    Mesh const& mesh,
    Field<real**> const& vf_field)
{
  std::stringstream stream;
  int const block = 0;
  dgt::vtk::write_vtr_start(stream, block, mesh, 0., 0);
  dgt::vtk::write_vtr_field(stream, "vfs", get_vfs(mesh, vf_field, block));
  dgt::vtk::write_vtr_end(stream);
  dgt::write_stream(path, stream);
}

static void insert_geom(
    mpicpp::comm* comm,
    std::string const& file,
    int const cells)
{
  auto const triangles_h = read(file);
  auto const triangles = to_device(triangles_h);
  auto const box = compute_bounding_box(triangles_h);
  Mesh mesh = make_mesh(comm, box, cells);
  auto const vfs = compute_vfs(mesh, triangles);
  write_mesh("debug.vtr", mesh, vfs);
}

int main(int argc, char** argv)
{
  dgt::initialize(argc, argv);
  if (argc != 3) {
    print_usage(argv[0]);
    return -1;
  }
  mpicpp::comm comm = mpicpp::comm::world();
  insert_geom(&comm, argv[1], std::stoi(argv[2]));
  dgt::finalize();
  return 0;
}
