#pragma once

#include <string>

#include "dgt_defines.hpp"
#include "dgt_grid.hpp"
#include "dgt_views.hpp"

namespace dgt {

struct FieldInfo {
  std::string name = "";
  int ent_dim = -1;
  int ncomps = -1;
};

class Field {
  private:
    FieldInfo m_info;
    View<double**> m_data[DIMS];
  public:
    Field() = default;
    void set_info(FieldInfo const& info);
    std::string name() const;
    int ncomps() const;
    int ent_dim() const;
    View<double**> data(int axis = 0) const;
    void allocate(p3a::grid3 const& cell_grid);
    void deallocate();
};

}
