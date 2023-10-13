#include <dgt_ghosting.hpp>
#include <dgt_partitioning.hpp>

#include <gtest/gtest.h>

using namespace dgt;
using namespace dgt::tree;

static Leaves refine_zleaf(
    int const dim,
    Leaves const& leaves,
    int const index)
{
  ZLeaves const z_leaves = order(dim, leaves);
  Marks marks(z_leaves.size(), REMAIN);
  marks[index] = REFINE;
  return modify(dim, z_leaves, marks);
}

static Leaves get_example_refined(int const dim)
{
  Leaves const leaves = create(dim, {3,3,0});
  return refine_zleaf(dim, leaves, 3);
}

static OwnedLeaves get_owned_leaves(
    mpicpp::comm* comm,
    ZLeaves const& z_leaves)
{
  using namespace linear_partitioning;
  int const rank = comm->rank();
  int const nranks = comm->size();
  int const nleaves = z_leaves.size();
  int const nlocal = get_num_local(nleaves, nranks, rank);
  int const offset = get_local_offset(nleaves, nranks, rank);
  auto begin = z_leaves.begin() + offset;
  auto end = z_leaves.begin() + offset + nlocal;
  return std::vector<ID>(begin, end);
}

TEST(ghosting, DEBUG)
{
  mpicpp::comm comm = mpicpp::comm::world();
  int const dim = 2;
  Vec3<bool> const periodic(false, false, false);
  Leaves const leaves = get_example_refined(dim);
  ZLeaves const z_leaves = order(dim, leaves);
  OwnedLeaves const owned_leaves = get_owned_leaves(&comm, z_leaves);
  Point const base_pt = get_base_point(dim, z_leaves);
  Adjacencies const adjs = get_adjacencies(dim, leaves, base_pt, periodic);
  Grid3 const cell_grid(6,6,0);
  ghosting::PackData data;
  data.build(cell_grid, adjs, owned_leaves, 5, 4);
}
