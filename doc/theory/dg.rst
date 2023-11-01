========
DG Bases
========

------------
Introduction
------------

Discontinuous Galerkin finite element basis functions can take many forms.
In particular, these basis functions are often categorized as *modal* or
*nodal*. Modal and nodal DG basis functions of the same polynomial order
typically span the same space, but differ in ways that we will not elaborate
on presently. In DGTile, we provide support for *modal* basis functions in
:math:`d=1,2,3` spatial dimensions of polynomial order
:math:`p=0,1,2,3` on a fixed reference cell :math:`\Omega^R := [-1,1]^d`.
These basis functions are constructed as products of Legendre polynomials
and are *orthogonal* over the reference cell.

----------------------
Spatial Representation
----------------------

Given a set of basis functions :math:`\{\phi_a\}_{a=1}^{n_b}`, the distribution
of a variable :math:`u` over an element :math:`\Omega^e` in the finite element
mesh can be expanded as:

.. math::
  u^e(\boldsymbol{\xi}, t) = \sum_{a=1}^{n_b} c_a^e(t) \phi_a(\boldsymbol{\xi}).

Here, :math:`\boldsymbol{\xi} \in [-1,1]^d` denotes a coordinate in the
reference cell, :math:`n_b` denotes the number of basis functions,
:math:`\phi_a` denotes an individual basis function that varies only in
space, and :math:`c_a^e(t)` represent *modal coefficients* that vary
only in time.

--------------------
Legendre Polynomials
--------------------

The first four Legendre polynomials of order :math:`p` on the one-dimensional
interval :math:`[-1,1]` are given as

.. math::
  \mathcal{L}_0(x) &= 1, \\
  \mathcal{L}_1(x) &= x, \\
  \mathcal{L}_2(x) &= \frac12 (3x^2-1), \\
  \mathcal{L}_3(x) &= \frac12 (5x^3-3x).

with derivatives

.. math::
  \mathcal{L}_0(x) &= 0, \\
  \mathcal{L}_1(x) &= 1, \\
  \mathcal{L}_2(x) &= 3x, \\
  \mathcal{L}_3(x) &= \frac32(5x^2-1).

The Legendre polynomials are *orthogonal* over the interval :math:`[-1,1]`,
where

.. math::
  \int_{-1}^1 \; \mathcal{L}_m(x) \mathcal{L}_n(x) \, \text{d}x =
  \frac{2}{2n + 1} \delta_{mn},

and :math:`\delta_{mn}` denotes the Kronecker delta.

---------------------
Construction of Bases
---------------------

In DGTile, we consider two families of DG basis functions, which
we refer to as *tensor* product bases and *complete polynomial*
bases, respectively.

Tensor Product Bases
--------------------

Tensor product bases are computed in two and
three dimensions as a tensor product of the Legendre polynomials,
meaning they include higher-order cross terms. In one, two, and
three spatial dimensions, the tensor product basis of order
:math:`p` can be written as:

.. math::
  \left\{ \phi_a(\boldsymbol{\xi}) \right\}_{a=1}^{n_b} &=
  \left\{ \mathcal{L}_i(\xi_1) \; | \; i \in \mathbb{N}_0 \land i \leq p
  \right\}, \\
  \left\{ \phi_a(\boldsymbol{\xi}) \right\}_{a=1}^{n_b} &=
  \left\{ \mathcal{L}_i(\xi_1) \mathcal{L}_j(\xi_2) \; | \; i,j \in
  \mathbb{N}_0 \land i,j \leq p \right\}, \\
  \left\{ \phi_a(\boldsymbol{\xi}) \right\}_{a=1}^{n_b} &=
  \left\{ \mathcal{L}_i(\xi_1) \mathcal{L}_j(\xi_2) \mathcal{L}_k(\xi_3) \; |
  \; i,j,k \in \mathbb{N}_0 \land i,j,k \leq p \right\},

respectively.

Complete Polynomial Bases
-------------------------

The complete polynomial bases are also constructed as products of
Legendre polynomials, but omit higher-order cross terms that are
included in a tensor product basis of the same order. In one,
two, and three spatial dimensions, the complete polynomial basis
of order :math:`p` can be written as:

.. math::
  \left\{ \phi_a(\boldsymbol{\xi}) \right\}_{a=1}^{n_b} &=
  \left\{ \mathcal{L}_i(\xi_1) \; | \; i \in \mathbb{N}_0 \land i \leq p
  \right\}, \\
  \left\{ \phi_a(\boldsymbol{\xi}) \right\}_{a=1}^{n_b} &=
  \left\{ \mathcal{L}_i(\xi_1) \mathcal{L}_j(\xi_2) \; | \; i,j \in
  \mathbb{N}_0 \land i + j \leq p \right\}, \\
  \left\{ \phi_a(\boldsymbol{\xi}) \right\}_{a=1}^{n_b} &=
  \left\{ \mathcal{L}_i(\xi_1) \mathcal{L}_j(\xi_2) \mathcal{L}_k(\xi_3) \; |
  \; i,j,k \in \mathbb{N}_0 \land i + j + k \leq p \right\},

respectively. Note that the one-dimensional tensor-product and complete
polynomial bases are identical.

Orthogonality
-------------

Simlar to the Legendre polynomials in one dimension, the tensor-product and
complete polynomial bases are *orthogonal* on the reference cell
:math:`\Omega^R`.

.. math::
  \int_{\Omega^R} \, \phi_a(\boldsymbol{\xi}) \phi_b(\boldsymbol{\xi}) \,
        \text{d} \Omega = \delta_{ab} M_a.

This is a desirable characteristic for explicit time-integration schemes,
as it makes the mass matrix :math:`M_a` *diagonal*, and thus readily
invertible.

.. note::

   TODO: We would like to document the values of the mass matrix here for the
   1/2/3D basis functions.
