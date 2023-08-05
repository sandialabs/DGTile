#include "dgt_when.hpp"

#include <gtest/gtest.h>

using namespace dgt;

static void test_when(
    WhenPtr when,
    std::vector<bool> const& expected_now,
    real const expected_final_time)
{
  real time = 0.;
  for (int step = 0; step < 6; ++step) {
    when->limit_time(time);
    bool const now = when->now(step, time);
    EXPECT_EQ(now, expected_now[step]);
    time += 0.1;
  }
  EXPECT_NEAR(time, expected_final_time, 1.e-15);
}

static WhenPtr combine_many_steps()
{
  WhenPtr when = std::make_shared<AtNever>();
  for (int step = 0; step < 6; ++step) {
    WhenPtr when2 = std::make_shared<AtStep>(step);
    when = combine_either(when, when2);
  }
  return when;
}

TEST(when, at_always)
{
  WhenPtr when = std::make_shared<AtAlways>();
  std::vector<bool> const expected = {1,1,1,1,1,1};
  test_when(when, expected, 0.6);
}

TEST(when, at_never)
{
  WhenPtr when = std::make_shared<AtNever>();
  std::vector<bool> const expected = {0,0,0,0,0,0};
  test_when(when, expected, 0.6);
}

TEST(when, at_step)
{
  WhenPtr when = std::make_shared<AtStep>(3);
  std::vector<bool> const expected = {0,0,0,1,0,0};
  test_when(when, expected, 0.6);
}

TEST(when, at_either_steps)
{
  WhenPtr a = std::make_shared<AtStep>(2);
  WhenPtr b = std::make_shared<AtStep>(4);
  WhenPtr when = std::make_shared<AtEither>(a, b);
  std::vector<bool> const expected = {0,0,1,0,1,0};
  test_when(when, expected, 0.6);
}

TEST(when, combine_many_steps)
{
  WhenPtr when = combine_many_steps();
  std::vector<bool> const expected = {1,1,1,1,1,1};
  test_when(when, expected, 0.6);
}

TEST(when, at_step_periodically)
{
  WhenPtr when = std::make_shared<AtStepPeriodically>(3);
  std::vector<bool> const expected = {1,0,0,1,0,0};
  test_when(when, expected, 0.6);
}

TEST(when, at_time) {
  WhenPtr when = std::make_shared<AtTime>(0.15);
  std::vector<bool> const expected = {0,0,1,0,0,0};
  test_when(when, expected, 0.6);
}

TEST(when, at_exact_time) {
  WhenPtr when = std::make_shared<AtExactTime>(0.15);
  std::vector<bool> const expected = {0,0,1,0,0,0};
  test_when(when, expected, 0.55);
}

TEST(when, at_time_periodically) {
  WhenPtr when = std::make_shared<AtTimePeriodically>(0.2);
  std::vector<bool> const expected = {1,0,1,0,1,0};
  test_when(when, expected, 0.6);
}

TEST(when, at_exact_time_periodically) {
  WhenPtr when = std::make_shared<AtExactTimePeriodically>(0.15);
  std::vector<bool> const expected = {1,0,1,0,1,0};
  test_when(when, expected, 0.5);
}
