#include "example.hpp"

namespace example {

static std::string get_banner()
{
  std::string banner;
  banner += " _____     ______     ______   __     __         ______    \n";
  banner += "/\\  __-.  /\\  ___\\   /\\__  _\\ /\\ \\   /\\ \\       /\\  ___\\   \n";
  banner += "\\ \\ \\/\\ \\ \\ \\ \\__ \\  \\/_/\\ \\/ \\ \\ \\  \\ \\ \\____  \\ \\  __\\   \n";
  banner += " \\ \\____-  \\ \\_____\\    \\ \\_\\  \\ \\_\\  \\ \\_____\\  \\ \\_____\\ \n";
  banner += "  \\/____/   \\/_____/     \\/_/   \\/_/   \\/_____/   \\/_____/ \n";
  return banner;
}

void run(State const& state)
{
  if (state.comm.rank() == 0) {
    printf("%s", get_banner().c_str());
    printf(" > running: '%s'\n", state.name.c_str());
    printf(" > from file: '%s'\n", state.input_file_name.c_str());
  }
}

}
