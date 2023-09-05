#include <dgt_initialize.hpp>

#include <gtest/gtest.h>

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  dgt::initialize(argc, argv);
  int const result = RUN_ALL_TESTS();
  dgt::finalize();
  return result;
}
