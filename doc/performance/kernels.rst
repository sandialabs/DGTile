====================
Kernel Consolidation
====================

------------
Introduction
------------

This page is intended to overview shared-memory parallelism design
decisions and their impacts on code performance on GPU architectures
only. As of this writing, the main branch of DGTile provides a
simple mechanism for shared-memory parallelism of physics kernels
in a tree-based AMR context. This simple method, however, incurs
a significant performance penalty when many AMR blocks are present
on a single GPU. A newer approach (currently in the higher-order
branch of FLEXO) addresses this performance penalty by launching
significantly fewer kernels. This new approach requires some C++
code sophistication and more advanced Kokkos usage techniques,
but seems worth the effort, as outlined in the results below.

---------------
Previous Method
---------------

In both approaches, a single GPU contains a collection of blocks
(axis-aligned cartesian grids) over which physics is performed.
In the previous approach for shared memory parallelism, a single
GPU kernel is launched *per block* every time its computation
is required. So, for instance, the computation of the numerical
fluxes in this approach launches `(3 * n_blocks)` kernels for
every time integration stage. Here the factor of `3` is
for the number of spatial dimensions. Additionally, for MPI
communication, a GPU kernel is launched for all 6 faces of each
block to perform the packing and unpacking of MPI message
information, resulting in `2 * 6 * n_blocks` GPU kernel launches
every time a transfer of border information is required.

----------
New Method
----------

In contrast, the new method provides shared-memory parallelism
over all blocks and all cells at the same time. This requires
two new abstractions that were not present in the previous
method. First, a dynamic data structure for the modal solution
information must be created, so that the entirety of the
solution information over each `n_block` blocks can be accessed
in a single kernel. Our chosen method to implement this dynamic
data structure is a `Kokkos::View` of `Kokkos::View`'s
(which we will endeavor to document more completely at a later
time). With this data structure in hand, we also require a
4D loop abstraction (in contrast to our previous 3D loop
abstractions), which loop over the indices i,j,k,b in the
Cartesian grid :math:`otimes` blocks product space. With these
in place, the corresponding numerical flux computation in
this approach launches `3` kernels for every time integration
stage (vs. `(3 * n_blocks)`). Additionally, we have
coalesced the packing and unpacking routines for MPI communication
into a single buffer, resulting in `2` GPU kernel launches
every time a transfer of border information is required.

-----------
Performance
-----------

As a very initial performance consideration, we compare timings
from the mini-app in the main branch of DGTile with timings
from the mini-app in the higher_order branch of DGTile for a
simple advection problem without any slope/positivity limiting
or visualization output. We consider an initial test case
of 4x4x4 blocks with 8x8x8 cells in each block on a single GPU
(determining an optimal block size will be an interesting topic
of future investigation). As a first comparison point, we will
simply consider the overall wall-time to solution for the two
approaches, as shown in the example outputs below:

.. code-block::
  :caption: Old approach

  ---
  dg hydro example!
  ---
   > name: compare
   > num mpi ranks: 1
   > polynomial order: 1
   > tensor product basis: 1
   > xmin: [0, 0, 0]
   > xmax: [1, 1, 1]
   > block grid: [4, 4, 4]
   > cell grid: [8, 8, 8]
   > periodic: [1, 1, 1]
   > init amr: 
   > initial conditions: advect
   > gamma: 1.4
   > final time: 1
   > CFL: 0.9
   > beta: 0.5
   > M: 1e+08
   > gravity: 0
   > gravity axis: 1
   > step frequency: 100
   > out frequency: 100000
   > amr frequency: -1
   > error regression: 1
  [step]       0 [t] 0.0000000e+00 [dt] 1.1310604e-03
  [step]     100 [t] 1.1309353e-01 [dt] 1.1309752e-03
  [step]     200 [t] 2.2618798e-01 [dt] 1.1309211e-03
  [step]     300 [t] 3.3928336e-01 [dt] 1.1308980e-03
  [step]     400 [t] 4.5237965e-01 [dt] 1.1309058e-03
  [step]     500 [t] 5.6547675e-01 [dt] 1.1309445e-03
  [step]     600 [t] 6.7857466e-01 [dt] 1.1310139e-03
  [step]     700 [t] 7.9167335e-01 [dt] 1.1311140e-03
  [step]     800 [t] 9.0477300e-01 [dt] 1.1310407e-03
  [step]     885 [t] 1.0000000e+00 [dt] 2.2253669e-04
  
  real  2m22.758s
  user  1m55.273s
  sys 0m26.731s

.. code-block::
  :caption: New approach

   _____     ______     ______   __     __         ______    
  /\  __-.  /\  ___\   /\__  _\ /\ \   /\ \       /\  ___\   
  \ \ \/\ \ \ \ \__ \  \/_/\ \/ \ \ \  \ \ \____  \ \  __\   
   \ \____-  \ \_____\    \ \_\  \ \_\  \ \_____\  \ \_____\ 
    \/____/   \/_____/     \/_/   \/_/   \/_____/   \/_____/ 
   > running: 'advection'
   > from file: 'test.lua'
  mesh stats:
  > blocks: 64
  > cells: 32768
  > minimum cell dx: [3.125000e-02, 3.125000e-02, 3.125000e-02]
  > maximum cell dx: [3.125000e-02, 3.125000e-02, 3.125000e-02]
  [step]: 0       [time]: 0.000000000000000e+00   [dt]: 1.131060390553168e-03
  [step]: 100     [time]: 1.130935266167213e-01   [dt]: 1.130975239244203e-03
  [step]: 200     [time]: 2.261879821031447e-01   [dt]: 1.130921098120344e-03
  [step]: 300     [time]: 3.392833622508727e-01   [dt]: 1.130897956940815e-03
  [step]: 400     [time]: 4.523796462452452e-01   [dt]: 1.130905767297175e-03
  [step]: 500     [time]: 5.654767549422777e-01   [dt]: 1.130944453457778e-03
  [step]: 600     [time]: 6.785746637451791e-01   [dt]: 1.131013912863365e-03
  [step]: 700     [time]: 7.916733491839930e-01   [dt]: 1.131114014544715e-03
  [step]: 800     [time]: 9.047729976757121e-01   [dt]: 1.131040717982266e-03
  
  real  0m40.100s
  user  0m27.506s
  sys 0m11.258s

That's over ``3.5`` times faster! Nice!

.. code-block::
  :caption: nsys profiling of old approach

  Time(%), Total Time (ns), Instances, Average,  Minimum, Maximum, Name
  24.4,    17075074573,     339840,    50244.5,  47456,   67712,   compute_intr_fluxes 
  23.0,    16112298547,     339840,    47411.4,  44032,   63327,   compute_side_integral 
  20.0,    14023823490,     113280,    123797.9, 122207,  158559,  compute_vol_integral 
  16.3,    11435717681,     679680,    16825.1,  15392,   27232,   fill_border 
  12.0,    8391473427,      679680,    12346.2,  11264,   21664,   compute_border_fluxes 
  3.0,     2125561930,      113280,    18763.8,  16799,   27616,   advance_explicitly 
  0.8,     587556176,       56640,     10373.5,  9311,    23616,   compute_dt 
  0.3,     199853143,       56640,     3528.5,   2880,    14464,   axpby 

.. code-block::
  :caption: nsys profiling of new approach

  Time(%), Total Time (ns), Instances, Average,   Minimum, Maximum, Name
  31.4,    10774213772,     5310,      2029042.1, 1746063, 2329737, compute_fluxes
  22.4,    7674387329,      1770,      4335812.1, 4281846, 4511636, pack 
  14.1,    4818854308,      1770,      2722516.6, 2503336, 2931747, volume_integral
  13.6,    4668608713,      1770,      2637632.0, 2539974, 2825796, face_integral
  6.6,     2271002011,      1770,      1283052.0, 1215444, 1351474, unpack
  5.3,     1812959779,      1770,      1024271.1, 890744,  1138965, zero_residual
  3.2,     1105485589,      1770,      624568.1,  551131,  709817,  advance_explicitly
  1.8,     619580668,       885,       700091.2,  636442,  800984,  axpby
  1.5,     517904319,       885,       585202.6,  467323,  670553,  compute_dt
