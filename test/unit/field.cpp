// a place where I am testing out some field stuff
// like the container abstraction for a modal field and how to
// parallelize over fields

#include <dgt_defines.hpp>
#include <dgt_macros.hpp>
#include <dgt_subgrid3.hpp>
#include <dgt_view.hpp>

#include <fmt/core.h>

#include <gtest/gtest.h>

using namespace dgt;

template <class Functor>
class Kokkos4DFunctor
{

  private:

    Functor m_functor;

  public:

    Kokkos4DFunctor(Functor functor) : m_functor(functor) {}

    DGT_ALWAYS_INLINE DGT_HOST_DEVICE inline
    auto operator()(
        int const b,
        int const i,
        int const j,
        int const k) const
    {
      return m_functor(b, Vec3<int>(i,j,k));
    }

};

template <class Functor>
void kokkos_4D_for_each(
    int const first,
    int const last,
    Vec3<int> const vec_first,
    Vec3<int> const vec_last,
    Functor functor)
{
  Vec3<int> const limits = vec_last - vec_first;
  if ((last - first) <= 0) return;
  if (limits.volume() == 0) return;
  using kokkos_policy = Kokkos::MDRangePolicy<
    Kokkos::IndexType<int>,
    Kokkos::Rank<4, Kokkos::Iterate::Left, Kokkos::Iterate::Left>>;
  Kokkos4DFunctor F(functor);
  Kokkos::parallel_for(
      "dgt::for_each(4D)",
      kokkos_policy(
        {vec_first.x(), vec_first.y(), vec_first.z(), first},
        {vec_last.x(), vec_last.y(), vec_last.z(), last}),
      F);
}

template <class Functor>
void for_each(
    int const num_blocks,
    Subgrid3 const& subgrid,
    Functor functor)
{
  kokkos_4D_for_each<Functor>(0, num_blocks, subgrid.lower(), subgrid.upper(), functor);
}

TEST(field, trying_stuff_out)
{

  int const num_blocks = 4;
  int const num_cells = 20;
  int const num_eqs = 10;
  int const num_modes = 8;

  std::vector<View<real***>> modal;
  modal.resize(num_blocks);
  for (int b = 0; b < num_blocks; ++b) {
    std::string const name = fmt::format("field_name(block={})", b);
    modal[b] = View<real***>(name, num_cells, num_eqs, num_modes);
  }

  for (int b = 0; b < num_blocks; ++b) {
    EXPECT_EQ(modal[b].extent(0), num_cells);
    EXPECT_EQ(modal[b].extent(1), num_eqs);
    EXPECT_EQ(modal[b].extent(2), num_modes);
  }

  HostPinnedView<View<real***>*> usable("test", num_blocks);
  for (int b = 0; b < num_blocks; ++b) {
    usable(b) = modal[0];
  }

  Subgrid3 cells({0,0,0}, {5,2,2});
  auto functor = [&] (int const block, Vec3<int> const& cell_ijk) {
    int const cell = cells.index(cell_ijk);
    usable[block](cell, 0, 0) = 1.;
  };
  for_each(num_blocks, cells, functor);

}
