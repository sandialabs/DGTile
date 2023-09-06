#include <dgt_lua.hpp>
#include <dgt_lua_interface.hpp>
#include <dgt_when.hpp>

#include <gtest/gtest.h>

using namespace dgt;

extern "C" {

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
    test_when(at_always, {1,1,1,1,1,1}, 0.6);
    test_when(at_never, {0,0,0,0,0,0}, 0.6);
    test_when(at_step, {0,0,0,1,0,0}, 0.6);
    test_when(at_either_steps, {0,0,1,0,1,0}, 0.6);
    test_when(at_many_steps, {1,1,1,1,1,1}, 0.6);
    test_when(at_step_periodically, {1,0,0,1,0,0}, 0.6);
    test_when(at_time, {0,0,1,0,0,0}, 0.6);
    test_when(at_exact_time, {0,0,1,0,0,0}, 0.55);
    test_when(at_time_periodically, {1,0,1,0,1,0}, 0.6);
    //test_when(at_exact_time_periodically, {1,0,1,0,1,0}, 0.5);
    return std::vector<lua::stack_object>({});
  });
}

}

TEST(when, from_lua)
{
  auto stack = lua::stack::newstate();
  stack.setglobal("dgt_check_whens", check_whens);
  std::filesystem::path const data_dir(std::getenv("INTEGRATION_DATA_DIR"));
  std::filesystem::path const when_path = data_dir / "when.lua";
  stack.dofile(when_path);
}
