#pragma once

#include <vector>

#include "p3a_grid3.hpp"
#include "p3a_simd_view.hpp"

#include "dgt_border.hpp"
#include "dgt_field.hpp"

namespace dgt {

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
    View<double***> m_path_cons[DIMS];
    View<double***> m_resid;
    std::vector<Field> m_fields;
  public:
    Block() = default;
    [[nodiscard]] int id() const;
    [[nodiscard]] int owner() const;
    [[nodiscard]] int dim() const;
    [[nodiscard]] int nsoln() const;
    [[nodiscard]] int nfields() const;
    [[nodiscard]] p3a::grid3 cell_grid() const;
    [[nodiscard]] p3a::box3<double> domain() const;
    [[nodiscard]] p3a::vector3<double> dx() const;
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
    [[nodiscard]] View<double***> path_cons(int axis) const;
    [[nodiscard]] View<double***> resid() const;
    [[nodiscard]] p3a::simd_view<double***> simd_soln(int idx) const;
    [[nodiscard]] p3a::simd_view<double***> simd_flux(int axis) const;
    [[nodiscard]] p3a::simd_view<double***> simd_path_cons(int axis) const;
    [[nodiscard]] p3a::simd_view<double***> simd_resid() const;
    [[nodiscard]] int field_idx(std::string name) const;
    [[nodiscard]] Field const& field(std::string name) const;
    void set_id(int id);
    void set_owner(int owner);
    void set_mesh(Mesh* mesh);
    void set_node(Node* node);
    void add_field(FieldInfo const& info);
    void reset();
    void allocate(int nsoln, int nmodal_eq, int nflux_eq);
    void deallocate();
};

}
