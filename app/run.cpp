#include <fstream>
#include <filesystem>

#include <dgt_print.hpp>

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
  real const css = total_cells * elapsed_steps / elapsed_time;
  previous_time = current_time;
  previous_step = current_step;
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

static void handle_step_output(Input const& in, State& state)
{
  auto const when = in.time.to_terminal;
  int const step = state.step;
  real const time = state.time;
  real const dt = state.dt;
  if (!when->now(step, time)) return;
  if (state.mesh.comm()->rank()) return;
  real const css = state.timer.compute(state.step);
  std::string const msg = fmt::format(
      "[step]: {:<7} [time]: {:.15e} [dt]: {:.15e} [cs/s]: {:.5e}",
      step, time, dt, css);
  printf("%s\n", msg.c_str());
}

static void handle_vtk_output(Input const& in, State& state)
{
  static int ctr = 0;
  auto const when = in.vtk.when;
  int const step = state.step;
  real const time = state.time;
  if (!when->now(step, time)) return;
  std::string const msg = fmt::format(
      "[visualization]: {:<7} [step]: {:<7} [time]: {:.15e}",
      ctr, step, time);
  printf("%s\n", msg.c_str());
  for (auto physics : state.physics) {
    physics->handle_visualization();
  }
  state.vtk_times.push_back(state.time);
  ctr++;
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

static void take_time_step(State& state)
{
  IntegratorPtr const& I = state.integrator;
  Physics& physics = state.physics;
  real const dt = state.dt;
  real const t = state.time;
  for (int stage = 0; stage < I->num_stages(); ++stage) {
    // TODO: might need integrator to handle the take time step
    // part so that different times / dts can be passed in here
    I->do_stage(physics, stage, dt, t);
  }
}

void run(mpicpp::comm* comm, Input const& in)
{
  print_banner(comm, in);
  create_output_dir(in.name);
  echo_input_file(in);
  State state;
  setup(comm, in, state);
  apply_initial_conditions(state);
  state.timer.update(state.mesh);
  while (true) {
    handle_vtk_output(in, state);
    if (state.time >= in.time.end_time) break;
    state.dt = compute_time_step(in, state);
    handle_step_output(in, state);
    take_time_step(state);
    state.time += state.dt;
    state.step++;
  }
  handle_step_output(in, state);
  write_pvd(in, state);
}

}
