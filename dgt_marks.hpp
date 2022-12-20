#pragma once

#include <vector>

namespace dgt {

class Mesh;

enum {REMAIN,REFINE,COARSEN};

std::vector<int8_t> reduce_marks(
    Mesh const& mesh,
    std::vector<int8_t> const& in_marks);

void ensure_2to1(Mesh const& mesh, std::vector<int8_t>& marks);

}
