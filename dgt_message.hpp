#pragma once

#include "mpicpp.hpp"

#include "dgt_views.hpp"

namespace dgt {

template <class T>
struct Message {
  public:
    View<T> val;
    mpicpp::request req;
  private:
    HostPinnedView<T> hval;
    bool needs_copy = false;
  public:
  void send(mpicpp::comm* comm, int rank, int tag) {
    needs_copy = false;
    copy(val, hval);
    req = comm->isend(hval.data(), hval.size(), rank, tag);
  }
  void recv(mpicpp::comm* comm, int rank, int tag) {
    needs_copy = true;
    resize(val, hval);
    req = comm->irecv(hval.data(), hval.size(), rank, tag);
  }
  void wait() {
    req.wait();
    if (needs_copy) {
      copy(hval, val);
      needs_copy = false;
    }
  }
};

}
