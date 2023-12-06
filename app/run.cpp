#include <fstream>
#include <filesystem>

#include "app.hpp"

namespace app {

static std::string get_banner()
{
  std::string banner;
  banner += " ┓\n";
  banner += "┏┫┏┓╋ ┏┓┏┓┏┓\n";
  banner += "┗┻┗┫┗ ┗┻┣┛┣┛\n";
  banner += "   ┛    ┛ ┛\n";
  return banner;
}

static void print_banner(mpicpp::comm* comm, Input const& in)
{
  if (comm->rank()) return;
  printf("%s", get_banner().c_str());
  printf("> running: '%s'\n", in.name.c_str());
  printf("> from file: '%s'\n", in.input_file_name.c_str());
}

static void create_output_dir(std::string const& name)
{
  std::filesystem::path const dir = name;
  std::filesystem::create_directory(dir);
}

void Timer::update(Mesh const& mesh)
{
  total_cells = mesh.num_total_cells();
}

real Timer::compute(int current_step)
{
  real const current_time = clock.seconds();
  real const elapsed_time = current_time - previous_time;
  real const elapsed_steps = current_step - previous_step;
  previous_time = current_time;
  previous_step = current_step;
  real const css = total_cells * elapsed_steps / elapsed_time;
  return css;
}

static void echo_input_file(Input const& in)
{
  std::filesystem::path const out_dir = in.name;
  std::filesystem::path const in_path = in.input_file_name;
  std::filesystem::path const echo_path = "input.echo.lua";
  std::filesystem::path const out_path = out_dir / echo_path;
  std::ifstream in_file(in_path);
  std::ofstream out_file(out_path);
  if (!in_file.is_open()) {
    throw std::runtime_error("-> couldn't open " + in_path.string());
  }
  if (!out_file.is_open()) {
    throw std::runtime_error("-> couldn't open " + out_path.string());
  }
  out_file << in_file.rdbuf();
  in_file.close();
  out_file.close();
}

static void apply_initial_conditions(State& state)
{
  for (auto physics : state.physics) {
    physics->apply_initial_conditions();
  }
}

static double compute_time_step(Input const& in, State& state)
{
  real min_dt = std::numeric_limits<double>::max();
  for (auto physics : state.physics) {
    real dt = physics->compute_time_step();
    min_dt = std::min(min_dt, dt);
  }
  min_dt = std::min(min_dt, in.time.end_time - state.time);
  return min_dt;
}

void run(mpicpp::comm* comm, Input const& in)
{
  print_banner(comm, in);
  create_output_dir(in.name);
  echo_input_file(in);
  State state;
  setup(comm, in, state);
  apply_initial_conditions(state);
  state.dt = compute_time_step(in, state);
}

}
