#include "dg.hpp"

TEST(dg_dim1, p0_q1_tensor)
{
  test_basis(1,0,1,true);
}

TEST(dg_dim1, p1_q2_tensor)
{
  test_basis(1,1,2,true);
}

TEST(dg_dim1, p2_q3_tensor)
{
  test_basis(1,2,3,true);
}

TEST(dg_dim1, p3_q4_tensor)
{
  test_basis(1,3,4,true);
}
