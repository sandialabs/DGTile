#include "dgt_amr.hpp"
#include "dgt_marks.hpp"
#include "dgt_mesh.hpp"
#include "dgt_spatial.hpp"

namespace dgt {

static void verify_marks(Mesh const& mesh, std::vector<int8_t> const& marks) {
  if (marks.size() != mesh.leaves().size()) {
    throw std::runtime_error("marks- invalid marks");
  }
}

std::vector<int8_t> reduce_marks(
    Mesh const& mesh,
    std::vector<int8_t> const& in_marks) {
  verify_marks(mesh, in_marks);
  mpicpp::comm* comm = mesh.comm();
  std::vector<int8_t> marks = in_marks;
  MPI_Allreduce(
      in_marks.data(), marks.data(), marks.size(),
      MPI_BYTE, MPI_MAX, comm->get());
  return marks;
}

}
