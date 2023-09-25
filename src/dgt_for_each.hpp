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
DGT_ALWAYS_INLINE DGT_HOST_DEVICE inline constexpr void inner_for_each(
    Subgrid3 const& s,
    Functor const& functor)
{
  for (int k = s.lower().z(); k < s.upper().z(); ++k) {
    for (int j = s.lower().y(); j < s.upper().y(); ++j) {
      for (int i = s.lower().x(); i < s.upper().x(); ++i) {
        functor(Vec3<int>(i, j, k));
      }
    }
  }
};

template <class Functor>
class Kokkos3DFunctor
{
  private:
    Functor m_functor;
  public:
    Kokkos3DFunctor(Functor const f) : m_functor(f) {}
    DGT_ALWAYS_INLINE DGT_HOST_DEVICE inline
    auto operator()(int const i, int const j, int const k) const
    {
      return m_functor(Vec3<int>(i,j,k));
    }
};

template <class Functor>
class Kokkos4DFunctor
{
  private:
    Functor m_functor;
  public:
    Kokkos4DFunctor(Functor const f) : m_functor(f) {}
    DGT_ALWAYS_INLINE DGT_HOST_DEVICE inline
    auto operator()(int const i, int const j, int const k, int const b) const
    {
      return m_functor(b, Vec3<int>(i,j,k));
    }
};

auto inline kokkos_3D_policy(
    Vec3<int> const& first,
    Vec3<int> const& last)
{
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
  return Kokkos::MDRangePolicy<
    Kokkos::IndexType<int>,
    Kokkos::Rank<3, Kokkos::Iterate::Left, Kokkos::Iterate::Left>>(
        {first.x(), first.y(), first.z()},
        {last.x(), last.y(), last.z()});
#else
  return Kokkos::MDRangePolicy<
    Kokkos::IndexType<int>,
    Kokkos::Rank<3, Kokkos::Iterate::Left, Kokkos::Iterate::Left>>(
        {first.x(), first.y(), first.z()},
        {last.x(), last.y(), last.z()},
        {1, 1, 1});
#endif
}

auto inline kokkos_4D_policy(
    int const sfirst,
    int const slast,
    Vec3<int> const& vfirst,
    Vec3<int> const& vlast)
{
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
  return Kokkos::MDRangePolicy<
    Kokkos::IndexType<int>,
    Kokkos::Rank<4, Kokkos::Iterate::Left, Kokkos::Iterate::Left>>(
        {vfirst.x(), vfirst.y(), vfirst.z(), sfirst},
        {vlast.x(), vlast.y(), vlast.z(), slast});
#else
  return Kokkos::MDRangePolicy<
    Kokkos::IndexType<int>,
    Kokkos::Rank<4, Kokkos::Iterate::Left, Kokkos::Iterate::Left>>(
        {vfirst.x(), vfirst.y(), vfirst.z(), sfirst},
        {vlast.x(), vlast.y(), vlast.z(), slast},
        {1, 1, 1, 1});
#endif
}

template <class Functor>
void kokkos_3D_for_each(
    std::string const& name,
    Vec3<int> const& vec_first,
    Vec3<int> const& vec_last,
    Functor const functor)
{
  Vec3<int> const vec_limits = vec_last - vec_first;
  if (vec_limits.volume() <= 0) return;
  auto policy_impl = kokkos_3D_policy(vec_first, vec_last);
  Kokkos3DFunctor functor_impl(functor);
  Kokkos::parallel_for(name, policy_impl, functor_impl);
}

template <class Functor>
void kokkos_4D_for_each(
    std::string const& name,
    int const scalar_first,
    int const scalar_last,
    Vec3<int> const& vec_first,
    Vec3<int> const& vec_last,
    Functor const functor)
{
  int const scalar_limits = scalar_last - scalar_first;
  Vec3<int> const vec_limits = vec_last - vec_first;
  if (scalar_limits <= 0) return;
  if (vec_limits.volume() <= 0) return;
  auto policy_impl = kokkos_4D_policy(
      scalar_first, scalar_last, vec_first, vec_last);
  Kokkos4DFunctor functor_impl(functor);
  Kokkos::parallel_for(name, policy_impl, functor_impl);
}

template <class Functor>
void for_each(
    std::string const& name,
    Subgrid3 const& s,
    Functor const f)
{
  int const dim = infer_dimension(s.upper());
  Subgrid3 const gs = generalize(dim, s);
  kokkos_3D_for_each(name, gs.lower(), gs.upper(), f);
}

template <class Functor>
void for_each(
    std::string const& name,
    int const n,
    Subgrid3 const& s,
    Functor const f)
{
  int const dim = infer_dimension(s.upper());
  Subgrid3 const gs = generalize(dim, s);
  kokkos_4D_for_each(name, 0, n, gs.lower(), gs.upper(), f);
}

}
