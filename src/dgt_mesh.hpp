#pragma once

#include <mpicpp.hpp>

#include "dgt_block_info.hpp"
#include "dgt_dg.hpp"
#include "dgt_tree.hpp"
#include "dgt_view.hpp"

namespace dgt {

class Mesh
{

  private:

    mpicpp::comm* m_comm = nullptr;

    Box3<real> m_domain = {{0,0,0}, {0,0,0}};
    Grid3 m_cell_grid = {-1,-1,-1};
    Vec3<bool> m_periodic = {false, false, false};

    Basis<View> m_basis;

    tree::Leaves m_leaves;
    tree::ZLeaves m_zleaves;
    std::vector<tree::ID> m_owned_leaves;

    BlockInfo m_block_info;

  public:

    Mesh() = default;

    void set_comm(mpicpp::comm* comm);
    void set_domain(Box3<real> const& domain);
    void set_cell_grid(Grid3 const& cell_grid);
    void set_periodic(Vec3<bool> const& periodic);
    void set_basis(Basis<View> m_basis);

    void verify();
    void init(Grid3 const& block_grid);

};

}
