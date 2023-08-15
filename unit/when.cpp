#include <dgt_lua.hpp>
#include <dgt_lua_interface.hpp>
#include <dgt_when.hpp>

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

static void test_at_always(WhenPtr when)
{
  std::vector<bool> const expected = {1,1,1,1,1,1};
  test_when(when, expected, 0.6);
}

static void test_at_never(WhenPtr when)
{
  std::vector<bool> const expected = {0,0,0,0,0,0};
  test_when(when, expected, 0.6);
}

static void test_at_step(WhenPtr when)
{
  std::vector<bool> const expected = {0,0,0,1,0,0};
  test_when(when, expected, 0.6);
}

static void test_at_either_steps(WhenPtr when)
{
  std::vector<bool> const expected = {0,0,1,0,1,0};
  test_when(when, expected, 0.6);
}

static void test_at_many_steps(WhenPtr when)
{
  std::vector<bool> const expected = {1,1,1,1,1,1};
  test_when(when, expected, 0.6);
}

static void test_at_step_periodically(WhenPtr when)
{
  std::vector<bool> const expected = {1,0,0,1,0,0};
  test_when(when, expected, 0.6);
}

static void test_at_time(WhenPtr when)
{
  std::vector<bool> const expected = {0,0,1,0,0,0};
  test_when(when, expected, 0.6);
}

static void test_at_exact_time(WhenPtr when)
{
  std::vector<bool> const expected = {0,0,1,0,0,0};
  test_when(when, expected, 0.55);
}

static void test_at_time_periodically(WhenPtr when)
{
  std::vector<bool> const expected = {1,0,1,0,1,0};
  test_when(when, expected, 0.6);
}

static void test_at_exact_time_periodically(WhenPtr when)
{
  std::vector<bool> const expected = {1,0,1,0,1,0};
  test_when(when, expected, 0.5);
}

TEST(when, at_always)
{
  WhenPtr when = std::make_shared<AtAlways>();
  test_at_always(when);
}

TEST(when, at_never)
{
  WhenPtr when = std::make_shared<AtNever>();
  test_at_never(when);
}

TEST(when, at_step)
{
  WhenPtr when = std::make_shared<AtStep>(3);
  test_at_step(when);
}

TEST(when, at_either_steps)
{
  WhenPtr a = std::make_shared<AtStep>(2);
  WhenPtr b = std::make_shared<AtStep>(4);
  WhenPtr when = std::make_shared<AtEither>(a, b);
  test_at_either_steps(when);
}

TEST(when, at_many_steps)
{
  WhenPtr when = combine_many_steps();
  test_at_many_steps(when);
}

TEST(when, at_step_periodically)
{
  WhenPtr when = std::make_shared<AtStepPeriodically>(3);
  test_at_step_periodically(when);
}

TEST(when, at_time) {
  WhenPtr when = std::make_shared<AtTime>(0.15);
  test_at_time(when);
}

TEST(when, at_exact_time) {
  WhenPtr when = std::make_shared<AtExactTime>(0.15);
  test_at_exact_time(when);
}

TEST(when, at_time_periodically) {
  WhenPtr when = std::make_shared<AtTimePeriodically>(0.2);
  test_at_time_periodically(when);
}

TEST(when, at_exact_time_periodically) {
  WhenPtr when = std::make_shared<AtExactTimePeriodically>(0.15);
  test_at_exact_time_periodically(when);
}

extern "C" {

static int check_whens(lua_State* L)
{
  return lua::function_wrapper(L,
  [] (lua::stack stack) {
    auto input = lua::table(stack.function_argument(1, "dgt_check_whens"));
    WhenPtr at_always = make_when(input.get_table(1));
    WhenPtr at_never = make_when(input.get_table(2));
    WhenPtr at_step = make_when(input.get_table(3));
    WhenPtr at_either_steps = make_when(input.get_table(4));
    WhenPtr at_many_steps = make_when(input.get_table(5));
    WhenPtr at_step_periodically = make_when(input.get_table(6));
    WhenPtr at_time = make_when(input.get_table(7));
    WhenPtr at_exact_time = make_when(input.get_table(8));
    WhenPtr at_time_periodically = make_when(input.get_table(9));
    WhenPtr at_exact_time_periodically = make_when(input.get_table(10));
    test_at_always(at_always);
    test_at_never(at_never);
    test_at_step(at_step);
    test_at_either_steps(at_either_steps);
    test_at_many_steps(at_many_steps);
    test_at_step_periodically(at_step_periodically);
    test_at_time(at_time);
    test_at_exact_time(at_exact_time);
    test_at_time_periodically(at_time_periodically);
    test_at_exact_time_periodically(at_exact_time_periodically);
    return std::vector<lua::stack_object>({});
  });
}

}

TEST(when, from_lua)
{
  auto stack = lua::stack::newstate();
  stack.setglobal("dgt_check_whens", check_whens);
  std::filesystem::path const data_dir(std::getenv("DATA_DIR"));
  std::filesystem::path const when_path = data_dir / "when.lua";
  stack.dofile(when_path);
}
