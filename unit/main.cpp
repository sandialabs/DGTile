#include <gtest/gtest.h>

#include <dgt_initialize.hpp>

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  dgt::initialize(argc, argv);
  int const result = RUN_ALL_TESTS();
  dgt::finalize();
  return result;
}
