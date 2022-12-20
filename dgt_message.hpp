#pragma once

#include "mpicpp.hpp"

namespace dgt {

template <class Data>
struct Message {
  Data val;
  mpicpp::request req;
  void send(mpicpp::comm* comm, int rank, int tag) {
    req = comm->isend(val.data(), val.size(), rank, tag);
  }
  void recv(mpicpp::comm* comm, int rank, int tag) {
    req = comm->irecv(val.data(), val.size(), rank, tag);
  }
  void wait() {
    req.wait();
  }
};

}
