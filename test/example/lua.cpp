#include <dgt_lua.hpp>
#include <dgt_lua_interface.hpp>

#include "example.hpp"

namespace example {

State make_state(
    dgt::lua::table const& in,
    std::string const& file_name)
{
  State result;
  result.comm = mpicpp::comm::world();
  result.input_file_name = file_name;
  result.name = in.get_string("name");
  result.basis = dgt::make_basis(in.get_table("basis"));
  return result;
}

}

extern "C" {

static int lua_dgtile_run(lua_State* L) {
  return dgt::lua::function_wrapper(L,
    [](dgt::lua::stack s) {
      auto dgtile_table = dgt::lua::table(s.getglobal("dgtile"));
      auto input_table = dgt::lua::table(s.function_argument(1, "dgtile.run"));
      std::string const file_name = dgtile_table.get_string("file");
      auto state = example::make_state(input_table, file_name);
      run(state);
      std::vector<dgt::lua::stack_object> results;
      return results;
    });
}

}

namespace example {

void load_dgtile_lua_module(
    dgt::lua::stack& stack,
    std::filesystem::path const& path)
{
  auto dgtile_table = stack.table("dgtile");
  dgtile_table.set("run", lua_dgtile_run);
  dgtile_table.set("file", path.string().c_str());
  stack.setglobal("dgtile", dgtile_table);
}

void run_lua_file(std::string const& path)
{
  auto stack = dgt::lua::stack::newstate();
  load_dgtile_lua_module(stack, path);
  stack.dofile(path);
}

}
