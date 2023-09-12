#include <fmt/core.h>

#include "dgt_cartesian.hpp"
#include "dgt_lua.hpp"
#include "dgt_lua_interface.hpp"

namespace dgt {

static void check_string_key(
    std::string const& function,
    lua::table_iterator::value_type const& pair,
    lua::table const& table)
{
  if (pair.first.type() != LUA_TSTRING) {
    std::string const msg = fmt::format(
        "dgt:{}[{}]-> key must be a string", function, table.name());
    throw std::runtime_error(msg);
  }
}

static void check_valid_when_keywords(lua::table const& in)
{
  for (auto pair : in) {
    check_string_key("make_when", pair, in);
    std::string const key = lua::string(pair.first).value();
    if (!((key == "kind") || (key == "frequency") || (key == "step") || (key == "time"))) {
      std::string const msg = fmt::format(
          "dgt:make_when[{}]-> unknown key `{}`", in.name(), key);
      printf("%s\n", msg.c_str());
      dgt::errors++;
    }
  }
}

static WhenPtr make_single_when(lua::table const& in)
{
  WhenPtr result = nullptr;
  check_valid_when_keywords(in);
  std::string const kind = in.get_string("kind");
  if (kind == "at_always") {
    result = std::make_shared<AtAlways>();
  } else if (kind == "at_never") {
    result = std::make_shared<AtNever>();
  } else if (kind == "at_step") {
    int const step = in.get_integer("step");
    result = std::make_shared<AtStep>(step);
  } else if (kind == "at_time") {
    real const time = in.get_number("time");
    result = std::make_shared<AtTime>(time);
  } else if (kind == "at_exact_time") {
    real const time = in.get_number("time");
    result = std::make_shared<AtExactTime>(time);
  } else if (kind == "at_step_periodically") {
    int const frequency = in.get_integer("frequency");
    result = std::make_shared<AtStepPeriodically>(frequency);
  } else if (kind == "at_time_periodically") {
    real const frequency = in.get_number("frequency");
    result = std::make_shared<AtTimePeriodically>(frequency);
  } else if (kind == "at_exact_time_periodically") {
    real const frequency = in.get_number("frequency");
    result = std::make_shared<AtExactTimePeriodically>(frequency);
  } else {
    std::string const msg = fmt::format(
        "dgt:when-> unknown kind `{}`", kind);
    printf("%s\n", msg.c_str());
    dgt::errors++;
  }
  return result;
}

WhenPtr make_when(lua::table const& in)
{
  WhenPtr result = std::make_shared<AtNever>();
  in.check_type(LUA_TTABLE);
  if (in.size() < 1) {
    std::string const msg = "dgt:make_when -> invalid table size";
    printf("%s\n", msg.c_str());
    dgt::errors++;
  }
  for (lua_Integer i = 1; i <= in.size(); ++i) {
    lua::table const when_table = in.get_table(i);
    WhenPtr when = make_single_when(when_table);
    if (when) result = combine_either(result, when);
  }
  return result;
}

static void check_valid_basis_keywords(lua::table const& in)
{
  for (auto pair : in) {
    check_string_key("make_basis", pair, in);
    std::string const key = lua::string(pair.first).value();
    if (!((key == "dimension") ||
          (key == "polynomial_order") ||
          (key == "quadrature_rule") ||
          (key == "tensor_product"))) {
      std::string const msg = fmt::format(
          "dgt::make_basis[{}]-> unknown key `{}`", in.name(), key);
      printf("%s\n", msg.c_str());
      dgt::errors++;
    }
  }
}

static void check_valid_basis_inputs(
    lua::table const& in,
    int const dim,
    int const p,
    int const q)
{
  if ((dim < 1) || (dim > 3)) {
    std::string const msg = fmt::format(
        "dgt::make_basis[{}]-> invalid dimension `{}`", in.name(), dim);
    throw std::runtime_error(msg);
  }
  if ((p < 0) || (p > max_polynomial_order)) {
    std::string const msg = fmt::format(
        "dgt::make_basis[{}]-> invalid polynomial_order `{}`", in.name(), p);
    throw std::runtime_error(msg);
  }
  if ((q < 1) || (q > max_1D_quadrature_points)) {
    std::string const msg = fmt::format(
        "dgt::make_basis[{}]-> invalid quadrature_rule `{}`", in.name(), q);
    throw std::runtime_error(msg);
  }
}


Basis<View> make_basis(lua::table const& in)
{
  check_valid_basis_keywords(in);
  int const dim = in.get_integer("dimension");
  int const p = in.get_integer("polynomial_order");
  int const q = in.get_integer("quadrature_rule");
  bool const tensor = in.get_boolean("tensor_product");
  check_valid_basis_inputs(in, dim, p, q);
  return build_basis<View>(dim, p, q, tensor);
}

static void check_valid_mesh_keywords(lua::table const& in)
{
  for (auto pair : in) {
    check_string_key("make_mesh", pair, in);
    std::string const key = lua::string(pair.first).value();
    if (!((key == "X") || (key == "Y") || (key == "Z"))) {
      std::string const msg = fmt::format(
          "dgt::make_mesh[{}]-> unknown key `{}`", in.name(), key);
      printf("%s\n", msg.c_str());
      dgt::errors++;
    }
  }
}

static void check_valid_axis_keywords(lua::table const& in)
{
  for (auto pair : in) {
    check_string_key("make_mesh", pair, in);
    std::string const key = lua::string(pair.first).value();
    if (!((key == "num_blocks") ||
          (key == "num_cells") ||
          (key == "min") ||
          (key == "max") ||
          (key == "periodic"))) {
      std::string const msg = fmt::format(
          "dgt::make_mesh[{}]-> unknown key `{}`", in.name(), key);
      printf("%s\n", msg.c_str());
      dgt::errors++;
    }
  }
}

Mesh make_mesh(
    lua::table const& in,
    Basis<View> const& basis,
    mpicpp::comm* comm)
{
  Mesh mesh;
  check_valid_mesh_keywords(in);
  Grid3 cell_grid(0,0,0);
  Grid3 block_grid(0,0,0);
  Vec3<bool> periodic(false, false, false);
  Box3<real> domain({0,0,0}, {0,0,0});
  for (int axis = 0; axis < DIMENSIONS; ++axis) {
    auto const axis_name = get_axis_name(axis);
    auto const axis_table = in.get_or_table(axis_name.c_str());
    check_valid_axis_keywords(axis_table);
    if (axis_table.num_entries() == 0) continue;
    cell_grid.extents()[axis] = axis_table.get_integer("num_cells");
    block_grid.extents()[axis] = axis_table.get_integer("num_blocks");
    domain.lower()[axis] = axis_table.get_number("min");
    domain.upper()[axis] = axis_table.get_number("max");
    periodic[axis] = axis_table.get_or("periodic", false);
  }
  mesh.set_comm(comm);
  mesh.set_domain(domain);
  mesh.set_cell_grid(cell_grid);
  mesh.set_periodic(periodic);
  mesh.set_basis(basis);
  mesh.init(block_grid);
  return mesh;
}

}
