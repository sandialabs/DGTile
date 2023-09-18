#pragma once

#include <Kokkos_DualView.hpp>

namespace dgt {

class Mesh;

template <class T>
using VtkView = Kokkos::DualView<T**, Kokkos::LayoutRight>;

void write_vtr_start(
    std::stringstream& stream,
    int const block,
    Mesh const& mesh,
    real const time,
    int const step);


void write_vtr_end(std::stringstream& stream);

}
