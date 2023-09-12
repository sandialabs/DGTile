#include <dgt_initialize.hpp>

#include <filesystem>

#include <fmt/core.h>

#include "example.hpp"

static void print_usage(std::string const& exe)
{
  printf("usage: %s <input.lua>\n", exe.c_str());
}

int main(int argc, char** argv)
{
  dgt::initialize(argc, argv);
  if (argc != 2) {
    print_usage(argv[0]);
    return -1;
  }
  try {
    std::string const input = argv[1];
    if (std::filesystem::exists(input)) {
      example::run_lua_file(input);
    } else {
      std::runtime_error(fmt::format("file {} doesn't exist", input));
    }
  } catch (std::exception const& e) {
    printf("DGTile: caught exception->\n%s\n", e.what());
    return -1;
  }
  dgt::finalize();
  return 0;
}
