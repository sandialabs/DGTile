#pragma once

#include <Kokkos_Core.hpp>

#include "dgt_grid3.hpp"
#include "dgt_subgrid3.hpp"

namespace dgt {

template <class Functor>
DGT_ALWAYS_INLINE inline constexpr void seq_for_each(
    Subgrid3 const& subgrid,
    Functor const& functor)
{
  int const dim = infer_dimension(subgrid.upper());
  Subgrid3 const s = generalize(dim, subgrid);
  if (s.size() == 0) return;
  for (int k = s.lower().z(); k < s.upper().z(); ++k) {
    for (int j = s.lower().y(); j < s.upper().y(); ++j) {
      for (int i = s.lower().x(); i < s.upper().x(); ++i) {
        functor(Vec3<int>(i, j, k));
      }
    }
  }
}

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
void for_each(int const n, Subgrid3 const& s, Functor f)
{
  int const dim = infer_dimension(s.upper());
  Subgrid3 const gs = generalize(dim, s);
  kokkos_4D_for_each(0, n, gs.lower(), gs.upper(), f);
}

}
