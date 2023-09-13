#pragma once

#include <mpicpp.hpp>

#include "dgt_block_info.hpp"
#include "dgt_dg.hpp"
#include "dgt_field.hpp"
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
    tree::OwnedLeaves m_owned_leaves;

    BlockInfo m_block_info;

    std::vector<SolutionField> m_modal;

  public:

    Mesh() = default;

    void set_comm(mpicpp::comm* comm);
    void set_domain(Box3<real> const& domain);
    void set_cell_grid(Grid3 const& cell_grid);
    void set_periodic(Vec3<bool> const& periodic);
    void set_basis(Basis<View> m_basis);

    void verify();

    void init(Grid3 const& block_grid);

    int modal_index(std::string const& name);
    void add_modal(
        std::string const& name,
        int const num_stored,
        int const num_comps,
        bool const with_flux = false);

    mpicpp::comm* comm() { return m_comm; }
    Box3<real> domain() const { return m_domain; }
    Grid3 cell_grid() const { return m_cell_grid; }
    Vec3<bool> periodic() const { return m_periodic; }
    Basis<View> const& basis() const { return m_basis; }
    tree::Leaves const& leaves() const { return m_leaves; }
    tree::ZLeaves const& z_leaves() const { return m_zleaves; }
    tree::OwnedLeaves const& owned_leaves() const { return m_owned_leaves; }
    BlockInfo const& block_info() const { return m_block_info; }

    int num_total_blocks() const;
    int num_owned_blocks() const;
    int num_total_cells() const;

    Field<real***> get_solution(std::string const& name, int const soln_idx);
    Vec3<Field<real***>> get_flux(std::string const& name);
    Field<real***> get_residual(std::string const& name);

};

}
