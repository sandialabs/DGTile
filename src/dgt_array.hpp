#pragma once

#include "dgt_defines.hpp"
#include "dgt_macros.hpp"

namespace dgt {

template <class T, int N>
class Array
{

  private:

    T m_data[N] = {};

  public:

    using reference = T&;
    using const_reference = T const&;

    DGT_METHOD constexpr reference operator[](int const i)
    {
      return m_data[i];
    }

    DGT_METHOD constexpr const_reference operator[](int const i) const
    {
      return m_data[i];
    }

};

}
