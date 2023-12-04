#include "dg.hpp"

static constexpr int dim = 2;

TEST(dg_dim2, p0_q1_complete)
{
  test_basis(dim,0,1,true);
}

TEST(dg_dim2, p1_q2_complete)
{
  test_basis(dim,1,2,false);
}

TEST(dg_dim2, p2_q3_complete)
{
  test_basis(dim,2,3,false);
}

TEST(dg_dim2, p3_q4_complete)
{
  test_basis(dim,3,4,false);
}

TEST(dg_dim2, p0_q1_tensor)
{
  test_basis(dim,0,1,true);
}

TEST(dg_dim2, p1_q2_tensor)
{
  test_basis(dim,1,2,true);
}

TEST(dg_dim2, p2_q3_tensor)
{
  test_basis(dim,2,3,true);
}

TEST(dg_dim2, p3_q4_tensor)
{
  test_basis(dim,3,4,true);
}
