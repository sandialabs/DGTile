#include "dgt_initialize.hpp"

int main(int argc, char** argv)
{
  dgt::initialize(argc, argv);
  dgt::finalize();
}
