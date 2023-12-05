#include <dgt_cartesian.hpp>
#include <dgt_lua.hpp>

#include "app.hpp"

namespace app {

static int parsing_errors = 0;

static void check_key_is_string(
    lua::table_iterator::value_type const& pair,
    lua::table const& in)
{
  if (pair.first.type() != LUA_TSTRING) {
    throw std::runtime_error(fmt::format(
          "{}-> key must be a string", in.name()));
  }
}

static void check_valid_keys(
    lua::table const& in,
    std::vector<std::string> const& keys)
{
  for (auto pair : in) {
    check_key_is_string(pair, in);
    std::string const key = lua::string(pair.first).value();
    if (std::count(keys.begin(), keys.end(), key)) continue;
    auto const msg =
      fmt::format("{}-> unknown key `{}`", in.name(), key);
    fprintf(stderr, "%s\n", msg.c_str());
    parsing_errors++;
  }
}

static void cond_err(
    lua::table const& in,
    std::string const& condition)
{
  auto const msg = fmt::format("{}.{}", in.name(), condition);
  fprintf(stderr, "%s\n", msg.c_str());
  parsing_errors++;
}

static WhenPtr make_single_when(lua::table const& in)
{
  WhenPtr result = nullptr;
  check_valid_keys(in, {
      "kind",
      "frequency",
      "step",
      "time"});
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
    auto const msg = fmt::format("{}-> unknown kind `{}`", in.name(), kind);
    fprintf(stderr, "%s\n", msg.c_str());
    parsing_errors++;
  }
  return result;
}

WhenPtr make_when(lua::table const& in)
{
  WhenPtr result = std::make_shared<AtNever>();
  in.check_type(LUA_TTABLE);
  if (in.size() < 1) {
    auto const msg = fmt::format("{}-> invalid table size", in.name());
    fprintf(stderr, "%s\n", msg.c_str());
    parsing_errors++;
  }
  for (lua_Integer i = 1; i <= in.size(); ++i) {
    lua::table const when_table = in.get_table(i);
    WhenPtr when = make_single_when(when_table);
    if (when) result = combine_either(result, when);
  }
  return result;
}

static inputs::Time parse_time(lua::table const& in)
{
  check_valid_keys(in, {
      "cfl",
      "end_time",
      "integrator",
      "to_terminal"});
  inputs::Time result;
  result.cfl = in.get_number("cfl");
  result.end_time = in.get_number("end_time");
  result.to_terminal = make_when(in.get_or_table("to_terminal"));
  result.integrator = in.get_string("integrator");
  if (result.cfl < 0) cond_err(in, "cfl < 0");
  if (result.end_time < 0) cond_err(in, "end_time < 0");
  if (result.integrator == "") cond_err(in, "integrator unset");
  return result;
}

static inputs::Basis parse_basis(lua::table const& in)
{
  check_valid_keys(in, {
      "polynomial_order",
      "quadrature_rule",
      "tensor_product"});
  inputs::Basis result;
  int const mp = max_polynomial_order;
  int const mq = max_1D_quadrature_points;
  std::string const smp = std::to_string(mp);
  std::string const smq = std::to_string(mq);
  int const p = in.get_integer("polynomial_order");
  int const q = in.get_integer("quadrature_rule");
  bool const tensor = in.get_boolean("tensor_product");
  result.polynomial_order = p;
  result.quadrature_rule = q;
  result.tensor_product = tensor;
  if (p < 0) cond_err(in, "polynomial_order < 0");
  if (p > mp) cond_err(in, "polynomial_order > " + smp);
  if (q < 1) cond_err(in, "quadrature_rule < 1");
  if (q > mq) cond_err(in, "quadrature_rule > " + smq);
  return result;
}

static void parse_mesh_axis(
    inputs::Mesh& result,
    lua::table const& in,
    int const axis)
{
  check_valid_keys(in, {
      "num_blocks",
      "num_cells",
      "min",
      "max",
      "periodic"});
  int const nblocks = in.get_integer("num_blocks");
  int const ncells = in.get_integer("num_cells");
  int const min = in.get_number("min");
  int const max = in.get_number("max");
  result.block_grid.extents()[axis] = nblocks;
  result.cell_grid.extents()[axis] = ncells;
  result.domain.lower()[axis] = min;
  result.domain.upper()[axis] = max;
  result.periodic[axis] = in.get_or("periodic", false);
  if (nblocks < 0) cond_err(in, "num_blocks < 0");
  if (ncells < 0) cond_err(in, "num_cells < 0");
  if (min > max) cond_err(in, "min > max");
}

static inputs::Mesh parse_mesh(lua::table const& in)
{
  check_valid_keys(in, {"X", "Y", "Z"});
  inputs::Mesh result;
  for (int axis = 0; axis < DIMENSIONS; ++axis) {
    auto const axis_name = get_axis_name(axis);
    auto axis_table = in.get_or_table(axis_name.c_str());
    if (axis_table.num_entries() == 0) continue;
    parse_mesh_axis(result, axis_table, axis);
  }
  return result;
}

Input make_input(
    lua::table const& in,
    std::string const& file_name)
{
  check_valid_keys(in, {
      "name",
      "time",
      "basis",
      "mesh"});
  Input result;
  result.name = in.get_string("name");
  result.input_file_name = file_name;
  result.time = parse_time(in.get_table("time"));
  result.basis = parse_basis(in.get_table("basis"));
  result.mesh = parse_mesh(in.get_table("mesh"));
  if (parsing_errors > 0) {
    throw std::runtime_error("dgtapp-> encountered parsing errors");
  }
  return result;
}

}

extern "C" {

static int lua_dgtile_run(lua_State* L) {
  return dgt::lua::function_wrapper(L, [](dgt::lua::stack s) {
    mpicpp::comm comm = mpicpp::comm::world();
    auto dgtile_table = dgt::lua::table(s.getglobal("dgtile"));
    auto input_table = dgt::lua::table(s.function_argument(1, "dgtile.run"));
    std::string const file_name = dgtile_table.get_string("file");
    app::Input const input = app::make_input(input_table, file_name);
    app::run(&comm, input);
    std::vector<dgt::lua::stack_object> results;
    return results;
  });
}

}

namespace app {

static void load_dgtile_lua_module(
    lua::stack& stack,
    std::filesystem::path const& path)
{
  auto dgtile_table = stack.table("dgtile");
  dgtile_table.set("run", lua_dgtile_run);
  dgtile_table.set("file", path.string().c_str());
  stack.setglobal("dgtile", dgtile_table);
}

void run_lua_file(std::string const& path)
{
  auto stack = lua::stack::newstate();
  load_dgtile_lua_module(stack, path);
  stack.dofile(path);
}

}
