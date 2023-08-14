#include <algorithm>
#include <limits>

#include <spdlog/spdlog.h>

#include "dgt_lua.hpp"
#include "dgt_when.hpp"

namespace dgt {

When::~When()
{
}

bool AtAlways::now(int const, real const)
{
  return true;
}

bool AtNever::now(int const, real const)
{
  return false;
}

AtStep::AtStep(int const step) :
  m_step(step)
{
}

bool AtStep::now(int const step, real const)
{
  return (step == m_step);
}

AtStepPeriodically::AtStepPeriodically(int const frequency)
  :m_frequency(frequency)
  ,m_last_step(std::numeric_limits<int>::min())
{
}

bool AtStepPeriodically::now(int const step, real const)
{
  if (step >= (m_last_step + m_frequency)) {
    m_last_step = step;
    return true;
  }
  return false;
}

AtTime::AtTime(real const time)
  :m_time(time)
  ,m_has_occurred(false)
{
}

bool AtTime::now(int const, real const time)
{
  if (m_has_occurred) {
    return false;
  } else {
    if (time >= m_time) {
      m_has_occurred = true;
      return true;
    } else {
      return false;
    }
  }
}

AtExactTime::AtExactTime(real const time)
  :AtTime(time)
{
}

void AtExactTime::limit_time(real& new_time) const
{
  if (!m_has_occurred) {
    new_time = std::min(m_time, new_time);
  }
}

AtTimePeriodically::AtTimePeriodically(real const frequency)
  :m_frequency(frequency)
  ,m_iteration(0)
{
}

bool AtTimePeriodically::now(int const, real const time)
{
  real next_time = m_iteration * m_frequency;
  if (time >= next_time) {
    while (time >= next_time) {
      ++m_iteration;
      next_time = m_iteration * m_frequency;
    }
    return true;
  } else {
    return false;
  }
}

AtExactTimePeriodically::AtExactTimePeriodically(real const frequency)
  :AtTimePeriodically(frequency)
{
}

void AtExactTimePeriodically::limit_time(real& new_time) const
{
  new_time = std::min(new_time, m_iteration * m_frequency);
}

AtEither::AtEither(WhenPtr first, WhenPtr second)
  :m_first(first)
  ,m_second(second)
{
}

void AtEither::limit_time(real& new_time) const
{
  m_first->limit_time(new_time);
  m_second->limit_time(new_time);
}

bool AtEither::now(int const step, real const time)
{
  bool const now_first = m_first->now(step, time);
  bool const now_second = m_second->now(step, time);
  return (now_first || now_second);
}

WhenPtr combine_either(WhenPtr a, WhenPtr b)
{
  WhenPtr when = std::make_shared<AtEither>(a, b);
  return when;
}

static void check_key(
    lua::table_iterator::value_type const& pair,
    lua::table const& table) {
  if (pair.first.type() != LUA_TSTRING) {
    spdlog::error("dgt:make_when[{}]-> key must be a string", table.name());
    throw std::runtime_error("dgt:make_when");
  }
}

static void check_valid_keywords(lua::table const& in)
{
  for (auto pair : in) {
    check_key(pair, in);
    std::string const key = lua::string(pair.first).value();
    if ((key != "kind") ||
        (key != "frequency") ||
        (key != "step") ||
        (key != "time")) {
      spdlog::error("dgt:make_when[{}]-> unknown key `{}`", in.name(), key);
      dgt::errors++;
    }
  }
}

static WhenPtr make_single_when(lua::table const& in)
{
  WhenPtr result = nullptr;
  check_valid_keywords(in);
  std::string const kind = in.get_string("kind");
  if (kind == "at_always") {
    result = std::make_shared<AtAlways>();
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
    spdlog::error("dgt:make_when -> invalid lua table");
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
