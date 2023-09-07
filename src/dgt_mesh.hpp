#pragma once

#include <mpicpp.hpp>

#include "dgt_dg.hpp"
#include "dgt_field.hpp"
#include "dgt_grid3.hpp"
#include "dgt_tree.hpp"

namespace dgt {

class Mesh
{

  private:

    mpicpp::comm* m_comm;

    Grid3 m_cell_grid;
    Vec3<bool> m_periodic;
    Box3<real> m_domain;

    Basis<View> m_basis;

    tree::Leaves m_leaves;
    tree::ZLeaves m_zleaves;
    tree::Leaves  m_owned_leaves;

    std::vector<Field> m_fields[Field::KINDS];

  public:

    Mesh() = default;

    void set_comm(mpicpp::comm* comm);
    void set_cell_grid(Grid3 const& cell_grid);
    void set_periodic(Vec3<bool> const& periodic);
    void set_domain(Box3<real> const& domain);
    void set_basis(int const dim, int const p, int const q, bool const tensor);
    void set_tree(tree::Leaves const& leaves);

};

}
