#include <spdlog/spdlog.h>

#include "dgt_lua.hpp"
#include "dgt_lua_interface.hpp"

namespace dgt {

static void check_when_key(
    lua::table_iterator::value_type const& pair,
    lua::table const& table) {
  if (pair.first.type() != LUA_TSTRING) {
    spdlog::error("dgt:make_when[{}]-> key must be a string", table.name());
    throw std::runtime_error("dgt:make_when");
  }
}

static void check_valid_when_keywords(lua::table const& in)
{
  for (auto pair : in) {
    check_when_key(pair, in);
    std::string const key = lua::string(pair.first).value();
    if (!((key == "kind") || (key == "frequency") || (key == "step") || (key == "time"))) {
      spdlog::error("dgt:make_when[{}]-> unknown key `{}`", in.name(), key);
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
    spdlog::error("dgt:when-> unknown kind `{}`", kind);
    dgt::errors++;
  }
  return result;
}

WhenPtr make_when(lua::table const& in)
{
  WhenPtr result = std::make_shared<AtNever>();
  in.check_type(LUA_TTABLE);
  if (in.size() < 1) {
    spdlog::error("dgt:make_when -> invalid table size");
    dgt::errors++;
  }
  for (lua_Integer i = 1; i <= in.size(); ++i) {
    lua::table const when_table = in.get_table(i);
    WhenPtr when = make_single_when(when_table);
    if (when) result = combine_either(result, when);
  }
  return result;
}

}
