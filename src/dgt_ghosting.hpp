#pragma once

namespace dgt {

class Mesh;

class Ghosting
{

  private:

    int m_num_messages = 0;

  public:

    void build(Mesh const& mesh);

};

}
