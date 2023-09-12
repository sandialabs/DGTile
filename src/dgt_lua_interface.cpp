#include <fmt/core.h>

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

Basis<View> make_basis(lua::table const& in)
{
  check_valid_basis_keywords(in);
  int const dim = in.get_integer("dimension");
  int const p = in.get_integer("polynomial_order");
  int const q = in.get_integer("quadrature_rule");
  bool const tensor = in.get_boolean("tensor_product");
  return build_basis<View>(dim, p, q, tensor);
}

}
