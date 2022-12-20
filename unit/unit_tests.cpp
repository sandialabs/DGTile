#include "gtest/gtest.h"

#include "dgt_library.hpp"

int main(int argc, char** argv) {
  int result = 0;
  {
    dgt::Library library(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
  }
  return result;
}
