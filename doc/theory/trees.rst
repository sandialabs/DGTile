=====
Trees
=====

------------
Introduction
------------

A *tree* is an abstract hierarchical data type that connects *nodes*
via *edges*. Each node in a tree is connected to a single *parent node*
with the exception of the *root node*, which has no parent.
Additionally, each node may be connected to multiple *child nodes*.
In DGTile, we use a *binary tree*, a *quadtree*, and an
*octree* to enable mesh adaptivity in one, two, and three spatial
dimensions, respectively. In a binary tree, quadtree, and octree each
*internal node* is connected to exactly two, four, and eight *child nodes*,
respectively. In contrast, a *leaf node* is any node in the tree that contains
no children. The *level* of a node in the tree is its distance away from
the rood node, where the *distance* is the number of edges along the
shortest path between two nodes. The root node in a tree has a level of 0.

--------------------------------------
Uniquely Identifying Nodes in the Tree
--------------------------------------

For the purposes of identifying, storing, modifying, and computing adjacencies
of leaves in the tree, it is useful to assign unique identifiers to every
possible node in the tree at every possible level. These identifiers are
computed using simple integer arithmetic. At every potential level :math:`l` in
the octree, we can assign a level-specific identifier :math:`\text{id}(l,i,j,k)`
to a potential node in the tree given the node's ijk Cartesian grid location
at that level. Presently, we consider the following choice for the
level-specific identifier

.. math::
  \text{id}(l,i,j,k) = i + 2^l(j + k 2^l)

where identifiers are strided the quickest along the x-axis and slowest
along the z-axis. Note that for a two-dimensional quadtree k will always
be 0 and likewise, for a one-dimensional binary tree k and j will always
be 0. The figure illustrates level-specific identifiers
for levels :math:`l=0,1,2` for an octree, where the x-axis extends from the
left to right of the page, the y-axis extends from out to into the page, and the
z-axis extends from bottom to top of the page. This convention will be used
in subsequent figures as well. For compatability with C++, all indices
start from 0.
