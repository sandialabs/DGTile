.. _dependencies:

============
Dependencies
============

This page documents the required dependencies of `DGTile`, how to obtain
these dependencies, and examples of how to install these dependencies
from source.

---
fmt
---

`{fmt}` is a modern formatting library. It is used in `DGTile` to format string
messages.

The `{fmt}` repository is hosted at `<https://github.com/fmtlib/fmt>`_.

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

  #!/bin/bash

  cmake .. \
  -DFMT_TEST=OFF \
  -DFMT_DOC=OFF \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_CXX_STANDARD=17 \
  -DCMAKE_CXX_EXTENSIONS=OFF \
  -DCMAKE_INSTALL_PREFIX="../inst"

----------
googletest
----------

`GTest` is a test framework. It is used in `DGTile` to support unit testing
of all externally facing methods.

The `GTest` repository is hosted at `<https://github.com/google/googletest>`_.

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
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCMAKE_INSTALL_PREFIX="../inst"

------
Kokkos
------

`Kokkos` is a performance portability layer.
It is used in `DGTile` to aid in the development of
shared-memory parallel algorithms and deployment across
multiple HPC platforms.

The `Kokkos` repository is hosted at `<https://github.com/kokkos/kokkos>`_.

.. code-block::

  git clone git@github.com:kokkos/kokkos.git
  git checkout 1a3ea28

To build and install the code:

.. code-block::

  cd kokkos
  mdkir build
  cd build
  * edit config.sh *
  chmod +x config.sh
  ./config.sh
  make install

.. code-block::
  :caption: Example Serial CMake configuration script config.sh

  #!/bin/bash

  cmake .. \
    -DCMAKE_CXX_STANDARD=17 \
    -DKokkos_ENABLE_COMPILE_AS_CMAKE_LANGUAGE=ON \
    -DKokkos_ENABLE_CUDA=OFF \
    -DKokkos_ENABLE_CUDA_LAMBDA=OFF \
    -DKokkos_ENABLE_OPENMP=OFF \
    -DKokkos_ENABLE_SERIAL=ON \
    -DKokkos_ENABLE_HEADER_SELF_CONTAINMENT_TESTS=OFF \
    -DKokkos_ENABLE_DEPRECATED_CODE_3=OFF \
    -DCMAKE_CUDA_SEPARABLE_COMPILATION=OFF \
    -DCMAKE_INSTALL_PREFIX="../inst"

---
Lua
---

`Lua` is a small programming language written in C. It is used
in `DGTile` to provide input parsing and scripting for its end
applications.

The `Lua` repository we make use of is a mirror of the official `Lua`
repositor with an additional CMake build system, and is located
at `<https://gitlab.com/codelibre/lua/lua-cmake>`_.

To obtain the repository:

.. code-block::

  git clone git@gitlab.com/codelibre/lua/lua-cmake
  git checkout dc4896f

To build and install the code:

.. code-block::

  cd lua-cmake 
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
    -DLUA_LANGUAGE="CXX" \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCMAKE_INSTALL_PREFIX="../inst"

------
mpicpp
------

`mpicpp` is a simple C++ wrapper around the MPI C library standard.
It is used in `DGTile` to simplify MPI usage for distributed memory
parallelism.

The `mpicpp` repository is hosted at `<https://github.com/sandialabs/mpicpp>`_.

A working MPI C++ compiler wrapper will be required to install `mpicpp`,
which we assume will be available on the desired HPC platform. Below, we
assume that this compiler has been aliased to `mpicxx`.

To obtain the repository:

.. code-block::

  git clone git@github.com:sandialabs/mpicpp.git
  git checkout 7e33b6a

To build and install the code:

.. code-block::

  cd mpicpp
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
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCMAKE_CXX_COMPILER=mpicxx \
    -DCMAKE_INSTALL_PREFIX="../inst"

----
zlib
----

`zlib-ng` is a portable compression library written in C.
It is used in `DGTile` to compress vtk visualization data.

The `zlib-ng` repository is hosted at `<https://github.com/zlib-ng/zlib-ng>`_.

To obtain the repository:

.. code-block::

  git clone git@github.com:zlib-ng/zlib-ng.git
  git checkout 9fb955b

To build and install the code:

.. code-block::

  cd zlib-ng
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
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCMAKE_INSTALL_PREFIX="../inst"
