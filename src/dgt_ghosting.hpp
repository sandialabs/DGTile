#pragma once

#include "dgt_defines.hpp"
#include "dgt_message.hpp"

namespace dgt {

class Mesh;

class Ghosting
{

  private:

    std::vector<Message<real>> m_messages[DIRECTIONS];

  public:

    void build(Mesh const& mesh);

  private:

    void build_messages(Mesh const& mesh);

};

}
