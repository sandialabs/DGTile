#include <fstream>
#include <filesystem>

#include <dgt_print.hpp>

#include "example.hpp"

namespace example {

static std::string get_banner()
{
  std::string banner;
  banner += " _____     ______     ______   __     __         ______    \n";
  banner += "/\\  __-.  /\\  ___\\   /\\__  _\\ /\\ \\   /\\ \\       /\\  ___\\   \n";
  banner += "\\ \\ \\/\\ \\ \\ \\ \\__ \\  \\/_/\\ \\/ \\ \\ \\  \\ \\ \\____  \\ \\  __\\   \n";
  banner += " \\ \\____-  \\ \\_____\\    \\ \\_\\  \\ \\_\\  \\ \\_____\\  \\ \\_____\\ \n";
  banner += "  \\/____/   \\/_____/     \\/_/   \\/_/   \\/_____/   \\/_____/ \n";
  return banner;
}

static void create_output_dir(std::string const& name)
{
  std::filesystem::path const dir = name;
  std::filesystem::create_directory(dir);
}

static void echo_input_file(Input const& in)
{
  using fmt::format;
  std::filesystem::path const out_dir = in.name;
  std::filesystem::path const in_path = in.input_file_name;
  std::filesystem::path const echo_path = "input.echo.lua";
  std::filesystem::path const out_path = out_dir / echo_path;
  std::ifstream in_file(in_path);
  std::ofstream out_file(out_path);
  if (!in_file.is_open()) {
    throw std::runtime_error("example-> couldn't open" + in_path.string());
  }
  if (!out_file.is_open()) {
    throw std::runtime_error("example-> couldn't open" + out_path.string());
  }
  out_file << in_file.rdbuf();
  in_file.close();
  out_file.close();
}

static void handle_step_output(Input const& in, State const& state)
{
  if (state.mesh.comm()->rank()) return;
  auto const when = in.time.to_terminal;
  int const step = state.step;
  real const time = state.time;
  real const dt = state.dt;
  if (!when->now(step, time)) return;
  std::string const msg = fmt::format(
      "[step]: {:<7} [time]: {:.15e}   [dt]: {:.15e}", step, time, dt);
  printf("%s\n", msg.c_str());
}

static void handle_viz_output(Input const& in, State& state)
{
  auto const when = in.vtk.when;
  int const step = state.step;
  real const time = state.time;
  if (!when->now(step, time)) return;
  if (state.mesh.comm()->rank() == 0) {
    std::string const msg = fmt::format(
        "<visualiztion> [step]: {:<7} [time] {:.15e}", step, time);
    printf("%s\n", msg.c_str());
  }
  write_out(in, state, 0);
}

static void write_pvd(Input const& in, State const& state)
{
  std::filesystem::path const out_dir(in.name + "/vtk");
  std::filesystem::path const pvd_path(in.name + ".pvd");
  std::filesystem::path const out_path = out_dir / pvd_path;
  std::stringstream stream;
  stream << "<VTKFile type=\"Collection\" version=\"0.1\">\n<Collection>\n";
  for (size_t i = 0; i < state.vtk_times.size(); ++i) {
    std::string const vtm_path = std::to_string(i) + "/blocks.vtm";
    stream << "<DataSet timestep=\"" << state.vtk_times[i] << "\" part=\"0\" ";
    stream << "file=\"" << vtm_path << "\"/>\n";
  }
  stream << "</Collection>\n</VTKFile>\n";
  write_stream(out_path, stream);
}

void run(mpicpp::comm* comm, Input const& in)
{

  if (comm->rank() == 0) {
    printf("%s", get_banner().c_str());
    printf(" > running: '%s'\n", in.name.c_str());
    printf(" > from file: '%s'\n", in.input_file_name.c_str());
  }

  create_output_dir(in.name);
  echo_input_file(in);
  State state;
  setup(state, comm, in);

  while (true) {
    handle_step_output(in, state);
    handle_viz_output(in, state);
    if (state.time >= in.time.end_time) break;
    state.dt = compute_dt(in, state);

    // stage 0
    zero_residual(state);
    state.mesh.ghost("hydro", 0);
    compute_fluxes(state, 0);
    compute_volume_integral(state, 0);
    compute_face_integral(state);
    advance_explicitly(state, 0, 1, state.dt);

    // stage 1
    zero_residual(state);
    state.mesh.ghost("hydro", 1);
    compute_fluxes(state, 1);
    compute_volume_integral(state, 1);
    compute_face_integral(state);
    advance_explicitly(state, 1, 1, state.dt);

    // final convex combination
    auto& mesh = state.mesh;
    auto& r = mesh.get_solution("hydro", 0);
    auto const& x = mesh.get_solution("hydro", 0);
    auto const& y = mesh.get_solution("hydro", 1);
    axpby(state, r, 0.5, x, 0.5, y);

    // time update
    state.time += state.dt;
    state.step++;
  }
  write_pvd(in, state);

}

}
