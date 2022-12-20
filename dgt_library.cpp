#include "mpicpp.hpp"

#include "p3a_opts.hpp"

#include "Kokkos_Core.hpp"

#include "caliper/cali-manager.h"

#include "dgt_library.hpp"

namespace dgt {

using namespace mpicpp;

class Library::impl {
  private:
    environment m_mpi;
    Kokkos::ScopeGuard m_kokkos;
    cali::ConfigManager m_caliper;
  public:
    impl(int& argc, char**& argv) :
      m_mpi(&argc, &argv),
      m_kokkos(argc, argv)
    {
      p3a::opts opts;
      opts.add("caliper").expect_arguments(1);
      opts.parse(argc, argv, true);
      if (opts.has("caliper")) {
        m_caliper.add(opts.argument("caliper").c_str());
      }
      if (m_caliper.error()) {
        throw std::runtime_error("caliper error: " +
            m_caliper.error_msg());
      }
      m_caliper.start();
    }
    ~impl() {
      m_caliper.flush();
    }
};

Library::Library(int& argc, char**& argv) :
  pimpl(new impl(argc, argv)) {
}

Library::~Library() {
}

}
