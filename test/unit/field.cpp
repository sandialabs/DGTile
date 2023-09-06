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

  int const num_blocks = 4;
  int const num_cells = 20;
  int const num_eqs = 5;
  int const num_modes = 8;
  ModalField field("hydro", num_blocks, num_cells, num_eqs, num_modes);
  auto const f = field.get();

  auto functor = [=] DGT_DEVICE (int const block) {
    int const cell = 0;
    int const eq = 0;
    int const mode = 0;
    f[block](cell, eq, mode) = 1.;
  };
  Kokkos::parallel_for("for", num_blocks, functor);

}

TEST(field, field_test)
{
  do_field_test();
}
