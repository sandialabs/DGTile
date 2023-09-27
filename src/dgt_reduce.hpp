#pragma once

#include "dgt_for_each.hpp"

namespace dgt {

template <class T, class Functor>
class Kokkos4DRFunctor
{
  private:
    Functor m_functor;
  public:
    Kokkos4DRFunctor(Functor const f) : m_functor(f) {}
    DGT_ALWAYS_INLINE DGT_HOST_DEVICE inline
    auto operator()(int const i, int const j, int const k, int const b, T& result) const
    {
      return m_functor(b, Vec3<int>(i,j,k), result);
    }
};

template <class T, class Functor, class Op>
void kokkos_4D_reduce_for_each(
    std::string const& name,
    int const scalar_first,
    int const scalar_last,
    Vec3<int> const& vec_first,
    Vec3<int> const& vec_last,
    Functor const functor,
    Op const op)
{
  int const scalar_limits = scalar_last - scalar_first;
  Vec3<int> const vec_limits = vec_last - vec_first;
  if (scalar_limits <= 0) return;
  if (vec_limits.volume() <= 0) return;
  auto policy_impl = kokkos_4D_policy(
      scalar_first, scalar_last, vec_first, vec_last);
  Kokkos4DRFunctor<T, Functor> functor_impl(functor);
  Kokkos::parallel_reduce(name, policy_impl, functor_impl, op);
}

template <class T, class Functor, class Op>
void reduce_for_each(
    std::string const& name,
    int const n,
    Subgrid3 const& s,
    Functor const f,
    Op const op)
{
  int const dim = infer_dimension(s.upper());
  Subgrid3 const gs = generalize(dim, s);
  kokkos_4D_reduce_for_each<T, Functor, Op>(
      name, 0, n, gs.lower(), gs.upper(), f, op);
}

}
