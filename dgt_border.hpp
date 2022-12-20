#pragma once

#include <vector>

#include "p3a_static_array.hpp"

#include "dgt_defines.hpp"
#include "dgt_basis.hpp"
#include "dgt_message.hpp"
#include "dgt_views.hpp"

namespace dgt {

class Node;
class Mesh;

enum {STANDARD=0,COARSE_TO_FINE=1,FINE_TO_COARSE=2,BOUNDARY=3};

static constexpr int NBORDER_CHILD = num_child(DIMS-1);

struct AMRBorderData {
  Message<double***> child_soln[NBORDER_CHILD];
  Message<double**> child_avg_soln[NBORDER_CHILD];
  View<double****> soln;
  View<double***> avg_soln;
};

class Border {
  private:
    int m_axis = -1;
    int m_dir = -1;
    int m_type = -1;
    Node const* m_node = nullptr;
    Node const* m_adj = nullptr;
  private:
    Message<double***> m_soln[ndirs];
    Message<double**> m_avg_soln[ndirs];
    AMRBorderData m_amr[ndirs];
    View<double****> m_amr_flux;
  public:
    Border() = default;
    [[nodiscard]] int axis() const;
    [[nodiscard]] int dir() const;
    [[nodiscard]] int type() const;
    [[nodiscard]] Node const* node() const;
    [[nodiscard]] Node const* adj() const;
    [[nodiscard]] Message<double***>& soln(int msg_dir);
    [[nodiscard]] Message<double**>& avg_soln(int msg_dir);
    [[nodiscard]] AMRBorderData& amr(int msg_dir);
    [[nodiscard]] p3a::static_array<View<double***>, ndirs> soln() const;
    [[nodiscard]] p3a::static_array<View<double****>, ndirs> amr_soln() const;
    [[nodiscard]] View<double****> amr_flux() const;
    void set_axis(int axis);
    void set_dir(int dir);
    void set_type(int type);
    void set_node(Node* node);
    void set_adj(Node* adj);
    void reset();
    void allocate(int nmodal_eq, int nflux_eq);
    void deallocate();
};

void begin_border_transfer(Mesh& m, int soln_idx);
void end_border_transfer(Mesh& m);

}
