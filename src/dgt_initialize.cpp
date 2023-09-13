#include <Kokkos_Core.hpp>

#include <mpi.h>

#include "dgt_defines.hpp"
#include "dgt_initialize.hpp"

namespace dgt {

void initialize(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  Kokkos::InitializationSettings kokkos_opts;
  kokkos_opts.set_map_device_id_by("mpi_rank");
  Kokkos::initialize(kokkos_opts);
}

void finalize()
{
  Kokkos::finalize();
  MPI_Finalize();
}

}
