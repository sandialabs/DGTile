#include "gtest/gtest.h"

#include "dgt_library.hpp"

using namespace dgt;

int main(int argc, char** argv) {
  int result = 0;
  {
    Library library(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
  }
  return result;
}
