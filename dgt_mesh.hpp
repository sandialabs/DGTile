#pragma once

#include <vector>

#include "mpicpp.hpp"

#include "p3a_box3.hpp"

#include "dgt_basis.hpp"
#include "dgt_tree.hpp"

namespace dgt {

using namespace p3a;

class Mesh {
  private:
    mpicpp::comm* m_comm = nullptr;
    box3<double> m_domain = {{0.,0.,0.}, {0.,0.,0.}};
    vector3<bool> m_periodic = {false,false,false};
    grid3 m_cell_grid = {0,0,0};
    Basis m_basis;
    int m_neq = -1;
    int m_nsoln = -1;
    std::vector<Node*> m_leaves;
    std::vector<Node*> m_owned_leaves;
    Tree m_tree;
  public:
    Mesh() = default;
    [[nodiscard]] mpicpp::comm* comm() const;
    [[nodiscard]] Tree& tree();
    [[nodiscard]] Tree const& tree() const;
    [[nodiscard]] int dim() const;
    [[nodiscard]] box3<double> domain() const;
    [[nodiscard]] vector3<bool> periodic() const;
    [[nodiscard]] grid3 cell_grid() const;
    [[nodiscard]] Basis& basis();
    [[nodiscard]] Basis const& basis() const;
    [[nodiscard]] int nsoln() const;
    [[nodiscard]] int neq() const;
    [[nodiscard]] std::vector<Node*> const& leaves() const;
    [[nodiscard]] std::vector<Node*> const& owned_leaves() const;
    void set_comm(mpicpp::comm* comm);
    void set_domain(box3<double> const& domain);
    void set_periodic(vector3<bool> const& periodic);
    void set_cell_grid(grid3 const& cell_grid);
    void set_nsoln(int nsoln);
    void set_neq(int neq);
    void set_tree(Tree& tree);
    void init(grid3 const& block_grid, int p, bool tensor);
    void scale(double l);
    void rebuild();
    void verify();
    void allocate();
    void clean();
};

void init_leaves(
    Mesh* mesh,
    Tree& tree,
    vector3<bool> const& periodic,
    std::vector<Node*> const& leaves);

}
