#include <dgt_field.hpp>
#include <dgt_view.hpp>
#include <dgt_subgrid3.hpp>

#include <fmt/core.h>

#include <gtest/gtest.h>

using namespace dgt;

template <class Functor>
class Kokkos4DFunctor
{
  private:
    Functor m_functor;
  public:
    Kokkos4DFunctor(Functor f) : m_functor(f) {}
    DGT_ALWAYS_INLINE DGT_HOST_DEVICE inline
    auto operator()(int const i, int const j, int const k, int const b) const
    {
      return m_functor(b, Vec3<int>(i,j,k));
    }
};

template <class Functor>
void kokkos_4D_for_each(
    int const scalar_first,
    int const scalar_last,
    Vec3<int> const& vec_first,
    Vec3<int> const& vec_last,
    Functor functor)
{
  int const scalar_limits = scalar_last - scalar_first;
  Vec3<int> const vec_limits = vec_last - vec_first;
  if (scalar_limits <= 0) return;
  if (vec_limits.volume() <= 0) return;
  using kokkos_policy = Kokkos::MDRangePolicy<
    Kokkos::IndexType<int>,
    Kokkos::Rank<4, Kokkos::Iterate::Left, Kokkos::Iterate::Left>>;
  kokkos_policy policy_impl(
      {vec_first.x(), vec_first.y(), vec_first.z(), scalar_first},
      {vec_last.x(), vec_last.y(), vec_last.z(), scalar_last});
  Kokkos4DFunctor functor_impl(functor);
  Kokkos::parallel_for("dgt::for_each(4D)", policy_impl, functor_impl);
}

template <class Functor>
void for_each(
    int const n,
    Subgrid3 const& s,
    Functor f)
{
  //TODO: maybe dimensionalize / generalize here
  kokkos_4D_for_each(0, n, s.lower(), s.upper(), f);
}

void do_field_test()
{


  // sizes
  int const num_eqs = 5;
  int const num_modes = 8;
  int const num_blocks = 4;
  Grid3 const cell_grid(5,2,2);

  // create the data
  ModalField hydro("hydro", num_blocks, cell_grid.size(), num_eqs, num_modes);
  auto const field = hydro.get();

  // loop over the data and do some things
  auto functor = [=] DGT_DEVICE (int const block, Vec3<int> const& cell_ijk) {
    int const cell = cell_grid.index(cell_ijk);
    for (int eq = 0; eq < num_eqs; ++eq) {
      for (int mode = 0; mode < num_modes; ++mode) {
        field[block](cell, eq, mode) = 1.;
      }
    }
  };
  for_each(num_blocks, cell_grid, functor);


  // compute a sum of the data of a field over a block
  int result = 0;
  View<real***> view = field[2];
  real* data = view.data();
  auto functor2 = [=] DGT_DEVICE (int const i, int& sum) {
    sum += data[i];
  };
  Kokkos::parallel_reduce("sum", view.size(), functor2, result);

  std::cout << "RESULT: " << result  << "\n";
  std::cout << " > view size: " <<  view.size() << "\n";
  std::cout << " > view extent 0: " << view.extent(0) << "\n";
  std::cout << " > view extent 1: " << view.extent(1) << "\n";
  std::cout << " > view extent 2: " << view.extent(2) << "\n";

}

TEST(field, field_test)
{
  do_field_test();
}
