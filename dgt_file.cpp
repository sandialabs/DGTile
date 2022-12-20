#include <fstream>

#include "caliper/cali.h"

#include "dgt_file.hpp"

namespace dgt {

void write_stream(
    std::filesystem::path const& path,
    std::stringstream const& stream) {
  CALI_CXX_MARK_FUNCTION;
  std::ofstream file_stream(path.c_str());
  if (!file_stream.is_open()) {
    throw std::runtime_error(
        "write_stream - could not open: " + path.string());
  }
  file_stream << stream.rdbuf();
  file_stream.close();
}

}
