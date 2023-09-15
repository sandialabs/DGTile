#include <dgt_cartesian.hpp>
#include <dgt_lua.hpp>
#include <dgt_when.hpp>

#include "example.hpp"

namespace example {

static int parsing_errors = 0;

static void check_string_key(
    lua::table_iterator::value_type const& pair,
    lua::table const& in)
{
  if (pair.first.type() != LUA_TSTRING) {
    throw std::runtime_error(fmt::format(
          "example:[{}]-> key must be a string", in.name()));
  }
}

static void check_valid_keys(
    lua::table const& in,
    std::vector<std::string> const& keys)
{
  for (auto pair : in) {
    check_string_key(pair, in);
    std::string const key = lua::string(pair.first).value();
    if (std::count(keys.begin(), keys.end(), key)) continue;
    auto const msg =
      fmt::format("example[{}]-> unknown key `{}`", in.name(), key);
    printf("%s\n", msg.c_str());
    parsing_errors++;
  }
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
    auto const msg = fmt::format(
        "example[{}]-> unknown kind `{}`", in.name(), kind);
    printf("%s\n", msg.c_str());
    parsing_errors++;
  }
  return result;
}

WhenPtr make_when(lua::table const& in)
{
  WhenPtr result = std::make_shared<AtNever>();
  in.check_type(LUA_TTABLE);
  if (in.size() < 1) {
    auto const msg = fmt::format(
        "example[{}]-> invalid table size", in.name());
    printf("%s\n", msg.c_str());
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
      "to_terminal"});
  inputs::Time result;
  result.cfl = in.get_number("cfl");
  result.end_time = in.get_number("end_time");
  result.to_terminal = make_when(in.get_or_table("to_terminal"));
  return result;
}

static inputs::Basis parse_basis(lua::table const& in)
{
  check_valid_keys(in, {
      "polynomial_order",
      "quadrature_rule",
      "tensor_product"});
  inputs::Basis result;
  result.polynomial_order = in.get_integer("polynomial_order");
  result.quadrature_rule = in.get_integer("quadrature_rule");
  result.tensor_product = in.get_boolean("tensor_product");
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
  result.block_grid.extents()[axis] = in.get_integer("num_blocks");
  result.cell_grid.extents()[axis] = in.get_integer("num_cells");
  result.domain.lower()[axis] = in.get_number("min");
  result.domain.upper()[axis] = in.get_number("max");
  result.periodic[axis] = in.get_or("periodic", false);
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

struct lua_scalar : inputs::function<real>
{
  lua::function f;
  lua_scalar(lua::function const& f_in) : f(f_in) {}
  ~lua_scalar() override {}
  real operator()(Vec3<real> const& x) override {
    auto f_results = f(x.x(), x.y(), x.z());
    real const f_val = lua::number(std::move(f_results[0])).value();
    return f_val;
  }
};

struct lua_vector : inputs::function<Vec3<real>>
{
  lua::function f;
  lua_vector(lua::function const& f_in) : f(f_in) {}
  ~lua_vector() override {}
  Vec3<real> operator()(Vec3<real> const& x) override {
    auto f_results = f(x.x(), x.y(), x.z());
    real const f_x = lua::number(std::move(f_results[0])).value();
    real const f_y = lua::number(std::move(f_results[1])).value();
    real const f_z = lua::number(std::move(f_results[2])).value();
    return Vec3<real>(f_x, f_y, f_z);
  }
};

inputs::function_ptr<real> make_scalar_function(
    lua::table const& in,
    std::string const& key)
{
  auto object = in.get_optional(key.c_str());
  if (!object.has_value()) {
    auto const msg = fmt::format(
        "example[{}]-> key `{}` doesn't exist", in.name(), key);
    printf("%s\n", msg.c_str());
    parsing_errors++;
    return nullptr;
  }
  lua::stack_object lua_so = object.value();
  if (lua_so.type() == LUA_TFUNCTION) {
    auto f = lua::function(in.get(key.c_str()));
    return std::make_unique<lua_scalar>(f);
  } else {
    auto const msg = fmt::format(
        "example[{}]-> key `{}` isn't a lua function", in.name(), key);
    printf("%s\n", msg.c_str());
    parsing_errors++;
    return nullptr;
  }
}

inputs::function_ptr<Vec3<real>> make_vector_function(
    lua::table const& in,
    std::string const& key)
{
  auto object = in.get_optional(key.c_str());
  if (!object.has_value()) {
    auto const msg = fmt::format(
        "example[{}]-> key `{}` doesn't exist", in.name(), key);
    printf("%s\n", msg.c_str());
    parsing_errors++;
    return nullptr;
  }
  lua::stack_object lua_so = object.value();
  if (lua_so.type() == LUA_TFUNCTION) {
    auto f = lua::function(in.get(key.c_str()));
    return std::make_unique<lua_vector>(f);
  } else {
    auto const msg = fmt::format(
        "example[{}]-> key `{}` isn't a lua function", in.name(), key);
    printf("%s\n", msg.c_str());
    parsing_errors++;
    return nullptr;
  }
}

static inputs::InitialConditions parse_ics(
    lua::table const& in,
    int const nmats)
{
  std::vector<std::string> valid_keys = {"velocity"};
  for (int mat = 0; mat < nmats; ++mat) {
    valid_keys.push_back(fmt::format("density_{}", mat));
    valid_keys.push_back(fmt::format("pressure_{}", mat));
  }
  check_valid_keys(in, valid_keys);
  inputs::InitialConditions result;
  result.densities.resize(nmats);
  result.pressures.resize(nmats);
  result.velocity = make_vector_function(in, "velocity");
  for (int mat = 0; mat < nmats; ++mat) {
    auto const rho_name = fmt::format("density_{}", mat);
    auto const p_name = fmt::format("pressure_{}", mat);
    result.densities[mat] = make_scalar_function(in, rho_name);
    result.pressures[mat] = make_scalar_function(in, p_name);
  }
  return result;
}

inputs::Materials parse_materials(
    lua::table const& in,
    int const nmats)
{
  std::vector<std::string> valid_keys;
  for (int mat = 0; mat < nmats; ++mat) {
    valid_keys.push_back(fmt::format("gamma_{}", mat));
  }
  check_valid_keys(in, valid_keys);
  inputs::Materials result;
  result.gammas.resize(nmats);
  for (int mat = 0;  mat < nmats; ++mat) {
    auto const gamma_name = fmt::format("gamma_{}", mat);
    result.gammas[mat] = in.get_number(gamma_name.c_str());
  }
  return result;
}

Input make_input(
    lua::table const& in,
    std::string const& file_name)
{
  check_valid_keys(in, {
      "name",
      "num_materials",
      "time",
      "basis",
      "mesh",
      "materials",
      "initial_conditions"});
  Input result;
  result.input_file_name = file_name;
  result.name = in.get_string("name");
  result.num_materials = in.get_integer("num_materials");
  result.time = parse_time(in.get_table("time"));
  result.basis = parse_basis(in.get_table("basis"));
  result.mesh = parse_mesh(in.get_table("mesh"));
  result.materials = parse_materials(
      in.get_table("materials"), result.num_materials);
  result.ics = parse_ics(
      in.get_table("initial_conditions"), result.num_materials);
  if (parsing_errors > 0) {
    throw std::runtime_error("example-> encountered errors when parsing");
  }
  return result;
}

}

extern "C" {

static int lua_dgtile_run(lua_State* L) {
  return dgt::lua::function_wrapper(L,
    [](dgt::lua::stack s) {
      mpicpp::comm comm = mpicpp::comm::world();
      auto dgtile_table = dgt::lua::table(s.getglobal("dgtile"));
      auto input_table = dgt::lua::table(s.function_argument(1, "dgtile.run"));
      std::string const file_name = dgtile_table.get_string("file");
      auto state = example::make_input(input_table, file_name);
      run(&comm, state);
      std::vector<dgt::lua::stack_object> results;
      return results;
    });
}

}

namespace example {

void load_dgtile_lua_module(
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
