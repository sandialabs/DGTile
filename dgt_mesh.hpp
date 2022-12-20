#pragma once

#include <vector>

#include "mpicpp.hpp"

#include "p3a_box3.hpp"

#include "dgt_basis.hpp"
#include "dgt_field.hpp"
#include "dgt_tree.hpp"

namespace dgt {

class Mesh {
  private:
    mpicpp::comm* m_comm = nullptr;
    p3a::box3<double> m_domain = {{0.,0.,0.}, {0.,0.,0.}};
    p3a::vector3<bool> m_periodic = {false,false,false};
    p3a::grid3 m_cell_grid = {0,0,0};
    Basis m_basis;
    int m_nsoln = -1;
    int m_nmodal_eq = -1;
    int m_nflux_eq = -1;
    std::vector<Node*> m_leaves;
    std::vector<Node*> m_owned_leaves;
    std::vector<FieldInfo> m_fields;
    Tree m_tree;
  public:
    Mesh() = default;
    [[nodiscard]] mpicpp::comm* comm() const;
    [[nodiscard]] Tree& tree();
    [[nodiscard]] Tree const& tree() const;
    [[nodiscard]] int dim() const;
    [[nodiscard]] p3a::box3<double> domain() const;
    [[nodiscard]] p3a::vector3<bool> periodic() const;
    [[nodiscard]] p3a::grid3 cell_grid() const;
    [[nodiscard]] Basis& basis();
    [[nodiscard]] Basis const& basis() const;
    [[nodiscard]] int nsoln() const;
    [[nodiscard]] int nmodal_eq() const;
    [[nodiscard]] int nflux_eq() const;
    [[nodiscard]] std::vector<Node*> const& leaves() const;
    [[nodiscard]] std::vector<Node*> const& owned_leaves() const;
    [[nodiscard]] std::vector<FieldInfo> const& fields() const;
    void set_comm(mpicpp::comm* comm);
    void set_domain(p3a::box3<double> const& domain);
    void set_periodic(p3a::vector3<bool> const& periodic);
    void set_cell_grid(p3a::grid3 const& cell_grid);
    void set_nsoln(int nsoln);
    void set_nmodal_eq(int neq);
    void set_nflux_eq(int neq);
    void set_tree(Tree& tree);
    void add_field(std::string name, int ent_dim, int ncomps);
    void init(p3a::grid3 const& block_grid, int p, bool tensor);
    void scale(double l);
    void rebuild();
    void verify();
    void allocate();
    void clean();
};

void init_leaves(
    Mesh* mesh,
    Tree& tree,
    p3a::vector3<bool> const& periodic,
    std::vector<Node*> const& leaves);

}
