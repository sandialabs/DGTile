#pragma once

#include "dgt_field.hpp"
#include "dgt_tree.hpp"

namespace dgt {


class Mesh
{

  private:

    Vec3<bool> m_periodic;
    Box3<real> m_domain;

    Grid3 m_cell_grid;

    tree::Leaves m_leaves;
    tree::ZLeaves m_zleaves;
    tree::Leaves m_owned_leaves;

    std::vector<ModalField> m_modal_fields;

};

}
