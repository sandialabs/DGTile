#include <dgt_partitioning.hpp>

#include <gtest/gtest.h>

using namespace dgt;
using namespace dgt::tree;

TEST(partitioning, linear_get_num_local)
{
  using namespace dgt::linear_partitioning;
  EXPECT_EQ(get_num_local(3,2,0), 2);
  EXPECT_EQ(get_num_local(3,2,1), 1);
  for (int i = 0; i < 4; ++i) {
    EXPECT_EQ(get_num_local(16,4,i), 4);
  }
  EXPECT_EQ(get_num_local(9,4,0), 3);
  EXPECT_EQ(get_num_local(9,4,1), 2);
  EXPECT_EQ(get_num_local(9,4,2), 2);
  EXPECT_EQ(get_num_local(9,4,3), 2);
}

TEST(partitioning, linear_get_local_offset)
{
  using namespace dgt::linear_partitioning;
  EXPECT_EQ(get_local_offset(3,2,0), 0);
  EXPECT_EQ(get_local_offset(3,2,1), 2);
  EXPECT_EQ(get_local_offset(9,4,0), 0);
  EXPECT_EQ(get_local_offset(9,4,1), 3);
  EXPECT_EQ(get_local_offset(9,4,2), 5);
  EXPECT_EQ(get_local_offset(9,4,3), 7);
}

TEST(partitioning, check_empty_ranks_will_work)
{
  using namespace dgt::linear_partitioning;
  int const num_total = 2;
  int const num_parts = 4;
  std::vector<int> ids{1,2};
  EXPECT_EQ(get_local_offset(num_total, num_parts, 0), 0);
  EXPECT_EQ(get_num_local(num_total, num_parts, 0), 1);
  EXPECT_EQ(get_local_offset(num_total, num_parts, 1), 1);
  EXPECT_EQ(get_num_local(num_total, num_parts, 1), 1);
  int const part = 3;
  int const offset = get_local_offset(num_total, num_parts, part);
  int const num_local = get_num_local(num_total, num_parts, part);
  auto begin = ids.begin() + offset;
  auto end = ids.begin() + offset + num_local;
  std::vector<int> owned_ids(begin, end);
  EXPECT_EQ(owned_ids.size(), 0);
}

TEST(partitioning, linear_get_owned_leaves)
{
  using namespace dgt::linear_partitioning;
  ZLeaves ids = {1,2,9,28,3,5,17,12};
  ZLeaves ids_part0 = {1,2,9,28};
  ZLeaves ids_part1 = {3,5,17,12};
  EXPECT_EQ(get_owned_leaves(0, 2, ids), ids_part0);
  EXPECT_EQ(get_owned_leaves(1, 2, ids), ids_part1);
}

TEST(partitioning, linear_get_part_info_serial)
{
  using namespace dgt::linear_partitioning;
  tree::GlobalToZ inv_z_leaves;
  inv_z_leaves[0] = 0;
  inv_z_leaves[1] = 1;
  inv_z_leaves[2] = 2;
  inv_z_leaves[3] = 3;
  for (int i = 0; i < 4; ++i) {
    PartInfo p = get_part_info(1, i, inv_z_leaves);
    EXPECT_EQ(p.rank, 0);
    EXPECT_EQ(p.block, i);
  }
}

TEST(partitioning, linear_get_part_info_2rank)
{
  using namespace dgt::linear_partitioning;
  tree::GlobalToZ inv_z_leaves;
  inv_z_leaves[0] = 0;
  inv_z_leaves[1] = 1;
  inv_z_leaves[2] = 2;
  inv_z_leaves[3] = 3;
  EXPECT_EQ(get_part_info(2, 0, inv_z_leaves).rank, 0);
  EXPECT_EQ(get_part_info(2, 0, inv_z_leaves).block, 0);
  EXPECT_EQ(get_part_info(2, 1, inv_z_leaves).rank, 0);
  EXPECT_EQ(get_part_info(2, 1, inv_z_leaves).block, 1);
  EXPECT_EQ(get_part_info(2, 2, inv_z_leaves).rank, 1);
  EXPECT_EQ(get_part_info(2, 2, inv_z_leaves).block, 0);
  EXPECT_EQ(get_part_info(2, 3, inv_z_leaves).rank, 1);
  EXPECT_EQ(get_part_info(2, 3, inv_z_leaves).block, 1);
}

TEST(partitioning, linear_get_part_info_3rank)
{
  using namespace dgt::linear_partitioning;
  tree::GlobalToZ inv_z_leaves;
  inv_z_leaves[0] = 0;
  inv_z_leaves[1] = 1;
  inv_z_leaves[2] = 2;
  inv_z_leaves[3] = 3;
  EXPECT_EQ(get_part_info(3, 0, inv_z_leaves).rank, 0);
  EXPECT_EQ(get_part_info(2, 0, inv_z_leaves).block, 0);
  EXPECT_EQ(get_part_info(3, 1, inv_z_leaves).rank, 0);
  EXPECT_EQ(get_part_info(3, 1, inv_z_leaves).block, 1);
  EXPECT_EQ(get_part_info(3, 2, inv_z_leaves).rank, 1);
  EXPECT_EQ(get_part_info(3, 2, inv_z_leaves).block, 0);
  EXPECT_EQ(get_part_info(3, 3, inv_z_leaves).rank, 2);
  EXPECT_EQ(get_part_info(3, 3, inv_z_leaves).block, 0);
}

TEST(partitioning, linear_get_part_info_4rank)
{
  using namespace dgt::linear_partitioning;
  tree::GlobalToZ inv_z_leaves;
  inv_z_leaves[0] = 0;
  inv_z_leaves[1] = 1;
  inv_z_leaves[2] = 2;
  inv_z_leaves[3] = 3;
  EXPECT_EQ(get_part_info(4, 0, inv_z_leaves).rank, 0);
  EXPECT_EQ(get_part_info(4, 0, inv_z_leaves).block, 0);
  EXPECT_EQ(get_part_info(4, 1, inv_z_leaves).rank, 1);
  EXPECT_EQ(get_part_info(4, 1, inv_z_leaves).block, 0);
  EXPECT_EQ(get_part_info(4, 2, inv_z_leaves).rank, 2);
  EXPECT_EQ(get_part_info(4, 2, inv_z_leaves).block, 0);
  EXPECT_EQ(get_part_info(4, 3, inv_z_leaves).rank, 3);
  EXPECT_EQ(get_part_info(4, 3, inv_z_leaves).block, 0);
}

TEST(partitioning, linear_get_part_info_5rank)
{
  using namespace dgt::linear_partitioning;
  tree::GlobalToZ inv_z_leaves;
  inv_z_leaves[0] = 0;
  inv_z_leaves[1] = 1;
  inv_z_leaves[2] = 2;
  inv_z_leaves[3] = 3;
  EXPECT_EQ(get_part_info(5, 0, inv_z_leaves).rank, 0);
  EXPECT_EQ(get_part_info(5, 0, inv_z_leaves).block, 0);
  EXPECT_EQ(get_part_info(5, 1, inv_z_leaves).rank, 1);
  EXPECT_EQ(get_part_info(5, 1, inv_z_leaves).block, 0);
  EXPECT_EQ(get_part_info(5, 2, inv_z_leaves).rank, 2);
  EXPECT_EQ(get_part_info(5, 2, inv_z_leaves).block, 0);
  EXPECT_EQ(get_part_info(5, 3, inv_z_leaves).rank, 3);
  EXPECT_EQ(get_part_info(5, 3, inv_z_leaves).block, 0);
}

TEST(partitioning, linear_get_extra_case)
{
  using namespace dgt::linear_partitioning;
  tree::GlobalToZ inv_z_leaves;
  inv_z_leaves[0] = 0;
  inv_z_leaves[1] = 1;
  inv_z_leaves[2] = 2;
  inv_z_leaves[3] = 3;
  inv_z_leaves[4] = 4;
  EXPECT_EQ(get_part_info(2, 0, inv_z_leaves).rank, 0);
  EXPECT_EQ(get_part_info(2, 0, inv_z_leaves).block, 0);
  EXPECT_EQ(get_part_info(2, 1, inv_z_leaves).rank, 0);
  EXPECT_EQ(get_part_info(2, 1, inv_z_leaves).block, 1);
  EXPECT_EQ(get_part_info(2, 2, inv_z_leaves).rank, 0);
  EXPECT_EQ(get_part_info(2, 2, inv_z_leaves).block, 2);
  EXPECT_EQ(get_part_info(2, 3, inv_z_leaves).rank, 1);
  EXPECT_EQ(get_part_info(2, 3, inv_z_leaves).block, 0);
  EXPECT_EQ(get_part_info(2, 4, inv_z_leaves).rank, 1);
  EXPECT_EQ(get_part_info(2, 4, inv_z_leaves).block, 1);
}
