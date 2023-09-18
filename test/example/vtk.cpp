#include <filesystem>

#include <dgt_print.hpp>
#include <dgt_vtk.hpp>

#include "example.hpp"

namespace example {

static void write_mesh(
    std::filesystem::path const& path,
    State const& state,
    int const soln_idx)
{
  Mesh const& mesh = state.mesh;
  int const nblocks = mesh.num_owned_blocks();
  for (int block = 0; block < nblocks; ++block) {
    std::stringstream stream;
    std::filesystem::path const block_path = path / fmt::format("{}.vtr", block);
    dgt::vtk::write_vtr_start(stream, block, mesh, state.time, state.step);
    dgt::vtk::write_vtr_end(stream);
    dgt::write_stream(block_path, stream);
  }
  if (state.mesh.comm()->rank() == 0) {
    std::filesystem::path const vtm_path = path / "blocks.vtm";
    std::stringstream stream;
    dgt::vtk::write_vtm(stream, "", mesh.num_total_blocks());
    dgt::write_stream(vtm_path, stream);
  }


  (void)soln_idx;
}

void write_out(Input const& in, State const& state, int soln_idx)
{
  static int ctr = 0;
  std::filesystem::path const out_dir = in.name;
  std::filesystem::path const vtk_dir = out_dir / "vtk";
  std::filesystem::path const path = vtk_dir / std::to_string(ctr);
  std::filesystem::create_directory(vtk_dir);
  std::filesystem::create_directory(path);
  write_mesh(path, state, soln_idx);
  ctr++;

  (void)soln_idx;
}

}
