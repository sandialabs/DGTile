#include <dgt_initialize.hpp>

#include <gtest/gtest.h>

#include <spdlog/spdlog.h>

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  spdlog::set_level(spdlog::level::debug);
  dgt::initialize(argc, argv);
  int const result = RUN_ALL_TESTS();
  dgt::finalize();
  return result;
}
