#pragma once

#include <memory>

namespace dgt {

class Library {
  private:
    class impl;
    std::unique_ptr<impl> pimpl;
  public:
    Library(int& argc, char**& argv);
    Library(Library&&) = delete;
    Library& operator=(Library&&) = delete;
    Library(Library const&) = delete;
    Library& operator=(Library const&) = delete;
    ~Library();
};

}
