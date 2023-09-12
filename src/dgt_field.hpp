#pragma once

#include <fmt/core.h>

#include "dgt_defines.hpp"
#include "dgt_view.hpp"

namespace dgt {

template <class T>
class Field {

  public:

    using view_t = View<T>;
    using uview_t = UnmanagedView<T>;
    using storage_t = std::vector<view_t>;
    using accessor_t = HostPinnedView<uview_t*>;

  protected:

    std::string m_name;
    storage_t m_storage;
    accessor_t m_accessor;

  public:

    Field() = default;

    std::string name() const { return m_name; }
    accessor_t get() { return m_accessor; }
    view_t get_view(int const block) { return m_storage[block]; }

    template <class Q = T>
    typename std::enable_if_t<View<Q>::rank == 1>
    create(
        std::string const& name,
        int const num_blocks,
        int const n0)
    {
      using Kokkos::view_alloc;
      using Kokkos::WithoutInitializing;
      m_name = name;
      m_storage.resize(num_blocks);
      m_accessor = accessor_t(name, num_blocks);
      for (int block = 0; block < num_blocks; ++block) {
        auto const block_name = fmt::format("{}[{}]", name, block);
        auto const alloc = view_alloc(block_name, WithoutInitializing);
        m_storage[block] = view_t(alloc, n0);
        m_accessor[block] = uview_t(m_storage[block].data(), n0);
      }
    }

    template <class Q = T>
    typename std::enable_if_t<View<Q>::rank == 2>
    create(
        std::string const& name,
        int const num_blocks,
        int const n0,
        int const n1)
    {
      using Kokkos::view_alloc;
      using Kokkos::WithoutInitializing;
      m_name = name;
      m_storage.resize(num_blocks);
      m_accessor = accessor_t(name, num_blocks);
      for (int block = 0; block < num_blocks; ++block) {
        auto const block_name = fmt::format("{}[{}]", name, block);
        auto const alloc = view_alloc(block_name, WithoutInitializing);
        m_storage[block] = view_t(alloc, n0, n1);
        m_accessor[block] = uview_t(m_storage[block].data(), n0, n1);
      }
    }

    template <class Q = T>
    typename std::enable_if_t<View<Q>::rank == 3>
    create(
        std::string const& name,
        int const num_blocks,
        int const n0,
        int const n1,
        int const n2)
    {
      using Kokkos::view_alloc;
      using Kokkos::WithoutInitializing;
      m_name = name;
      m_storage.resize(num_blocks);
      m_accessor = accessor_t(name, num_blocks);
      for (int block = 0; block < num_blocks; ++block) {
        auto const block_name = fmt::format("{}[{}]", name, block);
        auto const alloc = view_alloc(block_name, WithoutInitializing);
        m_storage[block] = view_t(alloc, n0, n1, n2);
        m_accessor[block] = uview_t(m_storage[block].data(), n0, n1, n2);
      }
    }

};

struct Modal
{
  std::string name = "";
  bool has_flux = true;
  Field<real***> modal;
  Field<real***> residual;
  Vec3<Field<real***>> flux;
};

}
