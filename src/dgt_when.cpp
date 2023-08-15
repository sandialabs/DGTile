#include <algorithm>
#include <limits>

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

}
