#pragma once

#include "mpicpp.hpp"

namespace dgt {

template <class T>
struct Message
{

  int rank = -1;
  int tag = -1;
  int size = -1;
  T* data = nullptr;
  mpicpp::request req;

  void send(mpicpp::comm* c) { req = c->isend(data, size, rank, tag); }
  void recv(mpicpp::comm* c) { req = c->irecv(data, size, rank, tag); }
  void wait() { req.wait(); }

};

}
