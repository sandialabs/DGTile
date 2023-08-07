#pragma once

#include <cfloat>

#include "dgt_defines.hpp"
#include "dgt_macros.hpp"

namespace dgt {

template <int N>
class Bounds
{

  private:

    real m_extrema[N][DIRECTIONS];

  public:

    DGT_HOST_DEVICE Bounds()
    {
      for (int i = 0; i < N; ++i) {
        m_extrema[i][MIN] = -DBL_MIN;
        m_extrema[i][MAX] =  DBL_MAX;
      }
    }

    DGT_METHOD inline real clamp(int const idx, real const in) const
    {
      if (in < m_extrema[idx][MIN]) return m_extrema[idx][MIN];
      if (in > m_extrema[idx][MAX]) return m_extrema[idx][MAX];
      return in;
    }

    DGT_METHOD inline real& operator()(int const idx, int const dir)
    {
      return m_extrema[idx][dir];
    }

    DGT_METHOD inline real const& operator()(int const idx, int const dir) const
    {
      return m_extrema[idx][dir];
    }

};

}
