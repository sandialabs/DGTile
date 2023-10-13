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
    Basis<HostView> m_basis_h;

    tree::Leaves m_leaves;
    tree::ZLeaves m_zleaves;
    tree::OwnedLeaves m_owned_leaves;
    tree::Adjacencies m_leaf_adjs;

    BlockInfo<View> m_block_info;
    BlockInfo<HostView> m_block_info_h;

    std::vector<ModalField> m_modal;

  public:

    Mesh() = default;

    void set_comm(mpicpp::comm* comm);
    void set_domain(Box3<real> const& domain);
    void set_cell_grid(Grid3 const& cell_grid);
    void set_periodic(Vec3<bool> const& periodic);
    void set_basis(int const p, int const q, bool const tensor);

    void init(Grid3 const& block_grid);

    void add_modal(
        std::string const& name,
        int const num_stored,
        int const num_comps,
        bool const with_flux = true);

    mpicpp::comm* comm() { return m_comm; }
    mpicpp::comm const* comm() const { return m_comm; }
    Box3<real> domain() const { return m_domain; }
    Grid3 cell_grid() const { return m_cell_grid; }
    Vec3<bool> periodic() const { return m_periodic; }
    tree::Leaves const& leaves() const { return m_leaves; }
    tree::ZLeaves const& z_leaves() const { return m_zleaves; }
    tree::OwnedLeaves const& owned_leaves() const { return m_owned_leaves; }
    tree::Adjacencies const& leaf_adjs() const { return m_leaf_adjs; }

    Basis<View> const& basis() const { return m_basis; }
    Basis<HostView> const& basis_h() const { return m_basis_h; }

    BlockInfo<View> const& block_info() const { return m_block_info; }
    BlockInfo<HostView> const& block_info_h() const { return m_block_info_h; }

    int dim() const;
    int num_total_blocks() const;
    int num_owned_blocks() const;
    int num_total_cells() const;
    int num_owned_cells() const;

    Field<real***>& get_solution(std::string const& name, int const soln_idx);
    Field<real***>& get_fluxes(std::string const& name, int const axis);
    Field<real***>& get_residual(std::string const& name);

    Field<real***> const& get_solution(std::string const& name, int const soln_idx) const;
    Field<real***> const& get_fluxes(std::string const& name, int const axis) const;
    Field<real***> const& get_residual(std::string const& name) const;

    void print_stats() const;

  private:

    void ensure_set();

};

}
