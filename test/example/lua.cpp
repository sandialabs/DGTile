#include <dgt_lua.hpp>

//static std::string get_banner()
//{
//  std::string banner;
//  banner += " _____     ______     ______   __     __         ______    \n";
//  banner += "/\\  __-.  /\\  ___\\   /\\__  _\\ /\\ \\   /\\ \\       /\\  ___\\   \n";
//  banner += "\\ \\ \\/\\ \\ \\ \\ \\__ \\  \\/_/\\ \\/ \\ \\ \\  \\ \\ \\____  \\ \\  __\\   \n";
//  banner += " \\ \\____-  \\ \\_____\\    \\ \\_\\  \\ \\_\\  \\ \\_____\\  \\ \\_____\\ \n";
//  banner += "  \\/____/   \\/_____/     \\/_/   \\/_/   \\/_____/   \\/_____/ \n";
//  return banner;
//}

extern "C" {

static int lua_dgtile_run(lua_State* L) {
  return dgt::lua::function_wrapper(L,
    [](dgt::lua::stack s) {
      auto dgtile_table = dgt::lua::table(s.getglobal("dgtile"));
      std::string const file_name = dgtile_table.get_string("file");
      (void)dgtile_table;
      (void)file_name;
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
