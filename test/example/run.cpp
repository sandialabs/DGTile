#include <fstream>
#include <filesystem>

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
      "[step]: {} [time]: {} [dt]: {}", step, time, dt);
  printf("%s\n", msg.c_str());
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

  // begin the loop here
  state.dt = compute_dt(in, state);
  handle_step_output(in, state);
  compute_fluxes(state, 0);
  compute_vol_integral(state, 0);

}

}
