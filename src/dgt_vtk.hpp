#pragma once

#include <Kokkos_DualView.hpp>

#include "dgt_defines.hpp"

namespace dgt {

class Mesh;

namespace vtk {

template <class T>
using VtkView = Kokkos::DualView<T**, Kokkos::LayoutRight>;

void write_vtr_start(
    std::stringstream& stream,
    int const block,
    Mesh const& mesh,
    real const time,
    int const step);

template <class T>
void write_vtr_field(
    std::stringstream& stream,
    std::string const& name,
    VtkView<T> f);

void write_vtr_end(std::stringstream& stream);

void write_vtm(
    std::stringstream& stream,
    std::string const& prefix,
    int const num_blocks);

}
}
