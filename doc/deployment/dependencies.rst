.. _dependencies:

============
Dependencies
============

This page documents the required dependencies of DGTile, how to obtain
these dependencies, and examples of how to install these dependencies
from source.

---
fmt
---

`{fmt}` is a modern formatting library. It is used in DGTile to format string
messages.

The `{fmt}` repository is hosted `here`_.

.. _here: https://github.com/fmtlib/fmt

To obtain the repository:

.. code-block::

  git clone git@github.com:fmtlib/fmt.git
  git checkout e150ea0

To build and install the code:

.. code-block::
  
  cd fmt
  mkdir build
  cd build
  * edit config.sh *
  chmod +x config.sh
  ./config.sh
  make install
  
.. code-block::
  :caption: Example CMake configuration script config.sh

  #!/bin/sh

  cmake .. \
  -DFMT_TEST=OFF \
  -DFMT_DOC=OFF \
  -DCMAKE_CXX_STANDARD=17 \
  -DCMAKE_CXX_EXTENSIONS=OFF \
  -DCMAKE_INSTALL_PREFIX="../inst"

----------
googletest
----------

`GTest` is a test framework. It is used in DGTile to support unit testing
of all externally facing methods.

The `GTest` repository is hosted `here`_.

.. _here: https://github.com/google/googletest

To obtain the repository:

.. code-block::

  git clone git@github.com:google/googletest.git
  git checkout 7e33b6a

To build and install the code:

.. code-block::

  cd googletest
  mdkir build
  cd build
  * edit config.sh *
  chmod +x config.sh
  ./config.sh
  make install

.. code-block::
  :caption: Example CMake configuration script config.sh

  #!/bin/bash

  cmake .. \
    -Dgtest_force_shared_crt=ON \
    -DCMAKE_INSTALL_PREFIX="../inst"
