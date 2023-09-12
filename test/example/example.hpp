#include <string>

#include <dgt_dg.hpp>
#include <dgt_view.hpp>

#include <mpicpp.hpp>

namespace example {

struct State
{
  mpicpp::comm comm;
  std::string name = "";
  std::string input_file_name = "";
  dgt::Basis<dgt::View> basis;
};

void run_lua_file(std::string const& path);
void run(State const& state);

}
