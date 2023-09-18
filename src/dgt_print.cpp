#include <fstream>
#include <stdexcept>

#include <fmt/core.h>

#include "dgt_print.hpp"

namespace dgt {

void write_stream(
    std::filesystem::path const& path,
    std::stringstream const& stream)
{
  std::ofstream file_stream(path.c_str());
  if (!file_stream.is_open()) {
    std::string const msg = fmt::format(
        "dgt::write_stream-> could not open {}\n", path.string());
    throw std::runtime_error(msg);
  }
  file_stream << stream.rdbuf();
  file_stream.close();
}

}
