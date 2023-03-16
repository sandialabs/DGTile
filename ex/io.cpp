#include <fstream>
#include <stdexcept>

#include "caliper/cali.h"

#include "dgt_file.hpp"
#include "dgt_grid.hpp"
#include "dgt_interp.hpp"
#include "dgt_print.hpp"

#include "hydro.hpp"

namespace hydro {

static void tokenize(
    std::string const& str,
    std::string const& delim,
    std::vector<std::string>& out) {
  out.resize(0);
  size_t start;
  size_t end = 0;
  while ((start = str.find_first_not_of(delim, end)) != std::string::npos) {
    end = str.find(delim, start);
    out.push_back(str.substr(start, end - start));
  }
}

static constexpr char const* WHITESPACE = " \n\r\t\f\v";

static std::string ltrim(const std::string& s) {
	size_t start = s.find_first_not_of(WHITESPACE);
	return (start == std::string::npos) ? "" : s.substr(start);
}

static std::string rtrim(const std::string& s) {
	size_t end = s.find_last_not_of(WHITESPACE);
	return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

static std::string trim(std::string const& s) {
	return rtrim(ltrim(s));
}

template <class T>
static p3a::vector3<T> parse_vec3(std::string const& s) {
  std::vector<std::string> ds;
  tokenize(s, ",", ds);
  p3a::vector3<T> val(T(0),T(0),T(0));
  val.x() = dgt::string_to_type<T>(trim(ds[0]));
  val.y() = dgt::string_to_type<T>(trim(ds[1]));
  val.z() = dgt::string_to_type<T>(trim(ds[2]));
  return val;
}

void parse_input(
    Input& in,
    mpicpp::comm* comm,
    std::string const& name) {
  in.comm = comm;
  in.name = name;
  std::ifstream file(name);
  if (!file.is_open()) {
    throw std::runtime_error("could not open: " + name);
  }
  std::string line;
  std::vector<std::string> dline;
  while (std::getline(file, line)) {
    if (line.rfind("#", 0) == 0) continue;
    if (std::all_of(line.begin(), line.end(), ::isspace)) continue;
    tokenize(line, ":", dline);
    assert(dline.size() == 2);
    std::string key = dline[0];
    std::string val = dline[1];
    key = trim(key);
    val = trim(val);
    if      (key == "p") in.p = dgt::string_to_type<int>(val);
    else if (key == "tensor") in.tensor = dgt::string_to_type<bool>(val);
    else if (key == "xmin") in.xmin = parse_vec3<double>(val);
    else if (key == "xmax") in.xmax = parse_vec3<double>(val);
    else if (key == "block_grid") in.block_grid = parse_vec3<int>(val);
    else if (key == "cell_grid") in.cell_grid = parse_vec3<int>(val);
    else if (key == "periodic") in.periodic = parse_vec3<bool>(val);
    else if (key == "init_amr") in.init_amr = val;
    else if (key == "amr") in.amr = val;
    else if (key == "ics") in.ics = val;
    else if (key == "gamma") in.gamma = dgt::string_to_type<double>(val);
    else if (key == "tfinal") in.tfinal = dgt::string_to_type<double>(val);
    else if (key == "CFL") in.CFL = dgt::string_to_type<double>(val);
    else if (key == "beta") in.beta = dgt::string_to_type<double>(val);
    else if (key == "M") in.M = dgt::string_to_type<double>(val);
    else if (key == "gravity") in.gravity = dgt::string_to_type<double>(val);
    else if (key == "gravity_axis") in.gravity_axis = dgt::string_to_type<int>(val);
    else if (key == "step_frequency") in.step_frequency = dgt::string_to_type<int>(val);
    else if (key == "out_frequency") in.out_frequency = dgt::string_to_type<double>(val);
    else if (key == "amr_frequency") in.amr_frequency = dgt::string_to_type<double>(val);
    else if (key == "error_regression") in.error_regression = dgt::string_to_type<double>(val);
    else {
      throw std::runtime_error("invalid input key: " + key);
    }
  }
}

void print_input(Input const& in) {
  if (in.comm->rank() != 0) return;
  std::cout << "---\n";
  std::cout << "dg hydro example!\n";
  std::cout << "---\n";
  std::cout << " > name: " << in.name << "\n";
  std::cout << " > num mpi ranks: " << in.comm->size() << "\n";
  std::cout << " > polynomial order: " << in.p << "\n";
  std::cout << " > tensor product basis: " << in.tensor << "\n";
  std::cout << " > xmin: " << in.xmin << "\n";
  std::cout << " > xmax: " << in.xmax << "\n";
  std::cout << " > block grid: " << in.block_grid << "\n";
  std::cout << " > cell grid: " << in.cell_grid << "\n";
  std::cout << " > periodic: " << in.periodic << "\n";
  std::cout << " > init amr: " << in.init_amr << "\n";
  std::cout << " > initial conditions: " << in.ics << "\n";
  std::cout << " > gamma: " << in.gamma << "\n";
  std::cout << " > final time: " << in.tfinal << "\n";
  std::cout << " > CFL: " << in.CFL << "\n";
  std::cout << " > beta: " << in.beta << "\n";
  std::cout << " > M: " << in.M << "\n";
  std::cout << " > gravity: " << in.gravity << "\n";
  std::cout << " > gravity axis: " << in.gravity_axis << "\n";
  std::cout << " > step frequency: " << in.step_frequency << "\n";
  std::cout << " > out frequency: " << in.out_frequency << "\n";
  std::cout << " > amr frequency: " << in.amr_frequency << "\n";
  std::cout << " > error regression: " << in.error_regression << "\n";
}

void verify_input(Input const& in) {
  if ((in.p < 0) || (in.p > dgt::max_p)) throw std::runtime_error("input - invalid p");
  if (in.gamma < 0) throw std::runtime_error("input - invalid gamma");
  if (in.tfinal <= 0) throw std::runtime_error("input - invalid final time");
  if ((in.CFL <= 0) || (in.CFL >= 1)) throw std::runtime_error("input - invalid CFL");
  if (in.step_frequency <= 0) throw std::runtime_error("input - invalid step frequency");
  if (in.out_frequency <= 0) throw std::runtime_error("input - invalid out frequency");
}

template <class T>
P3A_HOST_DEVICE void assign_variable(
    dgt::vtk::VizView<double> v, int idx, T const& val);

template <>
P3A_HOST_DEVICE void assign_variable(
    dgt::vtk::VizView<double> v, int idx, double const& val) {
  v.d_view(idx, 0) = val;
}

template <>
P3A_HOST_DEVICE void assign_variable(
    dgt::vtk::VizView<double> v, int idx, p3a::vector3<double> const& val) {
  v.d_view(idx, X) = val.x();
  v.d_view(idx, Y) = val.y();
  v.d_view(idx, Z) = val.z();
}

template <class Function>
dgt::vtk::VizView<double> get_variable(
    Function const& function,
    Block const& block,
    int num_comps,
    int soln_idx,
    double gamma = 0.) {
  Basis const b = block.basis();
  p3a::grid3 const g = block.cell_grid();
  p3a::grid3 const cell_grid = dgt::generalize(g);
  p3a::grid3 const inner_grid = dgt::tensor_bounds(b.dim, b.p);
  p3a::grid3 const viz_cell_grid = dgt::generalize(dgt::get_viz_cell_grid(g, b.p));
  View<double***> U = block.soln(soln_idx);
  dgt::vtk::VizView<double> var;
  Kokkos::resize(var, viz_cell_grid.size(), num_comps);
  auto f = [=] P3A_HOST_DEVICE (p3a::vector3<int> const& cell_ijk) {
    int const cell = cell_grid.index(cell_ijk);
    p3a::for_each(p3a::execution::hot, inner_grid,
      [&] (p3a::vector3<int> const& inner_ijk) P3A_ALWAYS_INLINE {
        int const pt = inner_grid.index(inner_ijk);
        p3a::vector3<int> const viz_cell_ijk = (b.p+1)*cell_ijk + inner_ijk;
        int const viz_cell = viz_cell_grid.index(viz_cell_ijk);
        auto const val = function(U, b, cell, pt, gamma);
        assign_variable(var, viz_cell, val);
      }
    );
  };
  p3a::for_each(p3a::execution::par, cell_grid, f);
  return var;
}

struct density {
  P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline double operator()(
      View<double***> U, dgt::Basis const& b, int cell, int pt, double) const {
    return interp_scalar_intr(U, b, cell, pt, RH);
  }
};

struct velocity {
  P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline p3a::vector3<double> operator()(
      View<double***> U, dgt::Basis const& b, int cell, int pt, double) const {
    double const rho = interp_scalar_intr(U, b, cell, pt, RH);
    p3a::vector3<double> const m = interp_vec3_intr(U, b, cell, pt, MM);
    return m/rho;
  }
};

struct pressure {
  P3A_ALWAYS_INLINE P3A_HOST_DEVICE inline double operator()(
      View<double***> U, dgt::Basis const& b, int cell, int pt, double gamma) const {
    p3a::static_vector<double, NEQ> const U_pt = dgt::interp_vec_intr<NEQ>(U, b, cell, pt);
    double const val = get_pressure(U_pt, gamma);
    return val;
  }
};

void write_mesh(
    std::filesystem::path const& path,
    State const& state,
    int soln_idx) {
  CALI_CXX_MARK_FUNCTION;
  Input const& in = state.in;
  Mesh const& mesh = state.mesh;
  std::filesystem::create_directory(path);
  density f_rho;
  velocity f_vel;
  pressure f_press;
  double const gamma = in.gamma;
  for (Node* leaf : mesh.owned_leaves()) {
    std::stringstream stream;
    Block const& block = leaf->block;
    std::filesystem::path block_path = path;
    block_path /= std::to_string(block.id()) + ".vtr";
    dgt::vtk::write_vtr_start(stream, block, state.t, state.step);
    dgt::vtk::write_field(stream, "density", get_variable(f_rho, block, 1, soln_idx));
    dgt::vtk::write_field(stream, "velocity", get_variable(f_vel, block, DIMS, soln_idx));
    dgt::vtk::write_field(stream, "pressure", get_variable(f_press, block, 1, soln_idx, gamma));
    dgt::vtk::write_vtr_end(stream);
    dgt::write_stream(block_path, stream);
  }
  if (state.mesh.comm()->rank() == 0) {
    std::filesystem::path vtm_path = path;
    vtm_path /= "blocks.vtm";
    std::stringstream stream;
    dgt::vtk::write_vtm(stream, "", mesh.leaves().size());
    dgt::write_stream(vtm_path, stream);
  }
}

void write_out(State& state, int soln_idx) {
  static int ctr = 0;
  double const freq = state.in.out_frequency;
  std::filesystem::path base(state.in.name);
  base = base.filename();
  base = std::filesystem::path(base.string() + "_viz");
  std::filesystem::create_directory(base);
  std::filesystem::path const path = base / ("out" + std::to_string(ctr));
  if (state.t >= freq * ctr) {
    write_mesh(path, state, soln_idx);
    state.out_times.push_back(state.t);
    ctr++;
  }
}

static void write_pvd_impl(
    State& state,
    std::string const& file_name,
    std::string const& parts_name) {
  std::filesystem::path base(state.in.name);
  base = base.filename();
  base = std::filesystem::path(base.string() + "_viz");
  std::filesystem::create_directory(base);
  std::filesystem::path const path = base / file_name;
  int const nout = state.out_times.size();
  std::stringstream stream;
  stream << "<VTKFile type=\"Collection\" version=\"0.1\">\n<Collection>\n";
  for (int i = 0; i < nout; ++i) {
    stream << "<DataSet timestep=\"" << state.out_times[i] << "\" part=\"0\" ";
    std::string vtm_path = "out" + std::to_string(i) + "/" + parts_name;
    stream << "file=\"" << vtm_path << "\"/>\n";
  }
  stream << "</Collection>\n</VTKFile>\n";
  dgt::write_stream(path, stream);
}

void write_pvd(State& state) {
  std::string const name = state.in.name + ".pvd";
  write_pvd_impl(state, name, "blocks.vtm");
}

}
