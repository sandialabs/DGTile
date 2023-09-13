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

void run(Input const& in)
{
  if (in.comm.rank() == 0) {
    printf("%s", get_banner().c_str());
    printf(" > running: '%s'\n", in.name.c_str());
    printf(" > from file: '%s'\n", in.input_file_name.c_str());
  }
}

}
