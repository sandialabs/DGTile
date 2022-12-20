#pragma once

#include <any>
#include <vector>

#include "p3a_grid3.hpp"
#include "p3a_simd_view.hpp"

#include "dgt_border.hpp"

namespace dgt {

using namespace p3a;

struct Basis;
class Mesh;
class Node;

class Block {
  private:
    int m_id = -1;
    int m_owner = -1;
    Mesh const* m_mesh = nullptr;
    Node const* m_node = nullptr;
  private:
    Border m_border[DIMS][ndirs];
    std::vector<View<double***>> m_soln;
    View<double***> m_flux[DIMS];
    View<double***> m_resid;
  public:
    std::any user_data;
  public:
    Block() = default;
    [[nodiscard]] int id() const;
    [[nodiscard]] int owner() const;
    [[nodiscard]] int dim() const;
    [[nodiscard]] int nsoln() const;
    [[nodiscard]] grid3 cell_grid() const;
    [[nodiscard]] box3<double> domain() const;
    [[nodiscard]] vector3<double> dx() const;
    [[nodiscard]] double cell_detJ() const;
    [[nodiscard]] double side_detJ(int axis) const;
    [[nodiscard]] double amr_side_detJ(int axis) const;
    [[nodiscard]] Mesh const* mesh() const;
    [[nodiscard]] Node const* node() const;
    [[nodiscard]] Basis const& basis() const;
    [[nodiscard]] Border& border(int axis, int dir);
    [[nodiscard]] Border const& border(int axis, int dir) const;
    [[nodiscard]] View<double***> soln(int idx) const;
    [[nodiscard]] View<double***> flux(int axis) const;
    [[nodiscard]] View<double***> resid() const;
    [[nodiscard]] simd_view<double***> simd_soln(int idx) const;
    [[nodiscard]] simd_view<double***> simd_flux(int axis) const;
    [[nodiscard]] simd_view<double***> simd_resid() const;
    void set_id(int id);
    void set_owner(int owner);
    void set_mesh(Mesh* mesh);
    void set_node(Node* node);
    void reset();
    void allocate(int nsoln, int neq);
    void deallocate();
};

}
