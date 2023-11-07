=======
Packing
=======

------------
Introduction
------------

Packing and unpacking MPI communication buffers is a critical component
of finite element codes and can have a great impact on the performancce
and scalability of application code. On this page, we document a simple
performance test included in the DGTile repository in the
`test/performance/packing.cpp` source file that mimics the types of
operations used to pack an MPI communication message using three
different approaches. Presently, we refer to these approaches as 
`Method A`, `Method B`, and `Method C`.

From a high level, this test performs the following steps. First, a large
global array is populated that is analagous to a double value stored on
every cell over every block on the current MPI rank. Then, a 'buffer' array
is initialized that corresponds to the data for every message that will be sent
by the current MPI rank. We then 'pack' this buffer by grabbing appropriate
values from the global cell array. In this simple profiling code, we assume
that every block on the MPI rank has a full set of 26 adjacent neighbors,
which corresponds to the use of a uniform grid in three spatial dimensions.
We are interested in seeing what the most appropriate packing method is for both
GPU and CPU performance for different problem sizes on a given MPI rank.

--------
Method A
--------

`Method A` describes possibly the simplest packing algorithm. Each block
is looped over and then each possible adjacency of an individual block is
looped over and an individual kernel to pack buffer data is then launched
in this inner loop.

.. code-block::
  :caption: Pseudocode for Method A

  for block : num blocks
    for adjacency : possible block adjacencies (26)
      subgrid adjacent_cells = grab_adjacent_cells_to_pack
      functor = (adjacenct_cell_id) {
        buffer_index = compute_buffer_index(adjacent_cell_id, block, adjacency)
        buffer[buffer_index] = global_cell_values(block, adjacent_cell_id)
      }
      parallel_for(adjacent_cells, functor)
    end
  end

--------
Method B
--------

`Method B` describes a variant of `Method A`, where we launch individual
kernels with so-called 'cuda streams' (when executing on GPUs).

.. code-block::
  :caption: Pseudocode for Method B

  weights = vec_of_num_adjacent_cells
  exec_spaces = partition_space_with_kokkos(weights)
  for block : num blocks
    for adjacency : possible block adjacencies (26)
      exec_idx = serialize(block, adjacency)
      subgrid adjacent_cells = grab_adjacent_cells_to_pack
      functor = (adjacenct_cell_id) {
        buffer_index = compute_buffer_index(adjacent_cell_id, block, adjacency)
        buffer[buffer_index] = global_cell_values(block, adjacent_cell_id)
      }
      parallel_for_with_cuda_streams(adjacent_cells, functor, exec_spaces[exec_idx])
    end
  end

--------
Method C
--------

`Method C` describes a packing alorithm where every cell on the current
MPI rank is looped over and then inside that loop, each individual cell
is checked to determine if its data should be packed or not. If not, the
for loop continues and if so, the cell data is placed in the appropriate
place in the MPI buffer array.

.. code-block::
  :caption: Pseudocode for Method C

  functor = (block_id, cell_id) {
    for adjacency:: possible block adjacencies (26)
      if ( don't need to pack cell_id ) continue
      buffer_index = compute_buffer_index(cell_id, block, adjacency)
      buffer[buffer_index] = global_cell_values(block, cell_id)
    end
  }
  parallel_for(all_blocks, all_cells_in_blocks, functor)

-------
Timings
-------

Here we compare timings obtained on my CPU laptop (2.3 GHz Quad-Core Intel Core i7)
for the various methods and on a V100 GPU machine. We performed the 'packing' operation
20 times in each test and measure the total time (in microseconds) for each method in
the tables below:

+----------+-----------------------+-----------------------+----------------------+----------------------+
| CPU      | 1 block [100,100,100] | 1 block [200,200,200] | 10 blocks [32,32,32] | 20 blocks [22,22,22] |
+==========+=======================+=======================+======================+======================+
| Method A | 21152                 | 75787                 | 9926                 | 10796                |
+----------+-----------------------+-----------------------+----------------------+----------------------+
| Method B | 13889                 | 70813                 | 8968                 | 8858                 |
+----------+-----------------------+-----------------------+----------------------+----------------------+
| Method C | 529073                | 3817890               | 196689               | 128529               |
+----------+-----------------------+-----------------------+----------------------+----------------------+

+----------+-----------------------+-----------------------+----------------------+----------------------+
| GPU      | 1 block [100,100,100] | 1 block [200,200,200] | 10 blocks [32,32,32] | 20 blocks [22,22,22] |
+==========+=======================+=======================+======================+======================+
| Method A | 2812                  | 4970                  | 23403                | 46849                |
+----------+-----------------------+-----------------------+----------------------+----------------------+
| Method B | 2900                  | 4480                  | 26247                | 52924                |
+----------+-----------------------+-----------------------+----------------------+----------------------+
| Method C | 2988                  | 18346                 | 1406                 | 1434                 |
+----------+-----------------------+-----------------------+----------------------+----------------------+
