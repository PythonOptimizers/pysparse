.. Introduction to Pysparse

========================
Introduction to Pysparse
========================

PySparse extends the Python interpreter by a set of sparse matrix types holding
double precision values. PySparse also includes modules that implement

- iterative Krylov methods for solving linear systems of equations,
- Diagonal (Jacobi) and SSOR preconditioners,
- interfaces to direct solvers for sparse linear systems of equations (SuperLU
  and UMFPACK),
- a Jacobi-Davidson eigenvalue solver for the symmetric, generalised matrix
  eigenvalue problem (JDSYM).

The above modules are implemented as C extension modules for maximum
performance.

PySparse uses `NumPy <http://numpy.scipy.org>`_ for handling dense vectors and
matrices and uses SuperLU and UMFPACK for factorizing general sparse matrices.

Module Overview
===============

``spmatrix``
------------

The ``spmatrix`` module is the foundation of the Pysparse package. It extends
the Python interpreter by three new types named ``ll_mat``, ``csr_mat`` and
``sss_mat``. These types represent sparse matrices in the LL-, the CSR- and
SSS-formats respectively (see :ref:`formats-page`). For all three formats,
double precision values (C type double) are used to represent the nonzero
entries.  The common way to use the spmatrix module is to first build a matrix
in the LL-format. The LL-matrix is manipulated until it has its final shape and
content. Afterwards it may be converted to either the CSR- or SSS-format, which
needs less memory and allows for fast matrix-vector multiplications. A ll_mat
object can be created from scratch, by reading data from a file (in
`MatrixMarket format <http://math.nist.gov/MatrixMarket>`_) or as a result of
matrix operation (as e.g. a matrix-matrix multiplication). The ll_mat object
supports manipulating (reading, writing, add-updating) single entries or
sub-matrices. On the other hand, ``csr_mat`` and ``sss_mat`` matrices are not
constructed directly, instead they are created by converting ``ll_mat``
objects. Once created, ``csr_mat`` and ``sss_mat`` objects cannot be
manipulated. Their purpose is to support efficient matrix-vector
multiplications.

``itsolvers``
-------------

The ``itsolvers`` module provides a set of iterative methods for solving linear
systems of equations. The iterative methods are callable like ordinary Python
functions. All these functions expect the same parameter list, and all function
return values also follow a common standard. Any user-defined iterative solvers
should also follow these conventions, since other PySparse modules rely on them
(e.g. the ``jdsym`` module).

Currently the ``itsolvers`` module contains the following iterative methods:
PCG, MINRES, QMRS, BICGSTAB and CGS.

``precon``
----------

The ``precon`` module provides preconditioners, which can be used e.g. for the
iterative methods implemented in the ``itsolvers`` module or the JDSYM
eigensolver (in the ``jdsym`` module).  In the PySparse framework, any Python object
that has the following properties can be used as a preconditioner:

- a ``shape`` attribute, which returns a 2-tuple describing the dimension of the
  preconditioner,
- a ``precon`` method, that accepts two vectors x and y, and applies the
  preconditioner to x and stores the result in y. Both x and y are double
  precision, rank-1 NumPy arrays of appropriate size.

The ``precon`` module currently implements m-step Jacobi and m-step SSOR
preconditioners.

``superlu``
-----------

The ``superlu`` module interfaces the `SuperLU
<http://crd.lbl.gov/~xiaoye/SuperLU/>`_ library to make it usable by Python
code. SuperLU is a software package written in C for the direct solution of
a general linear system of equations. SuperLU computes LU-factorizations of
general non-symmetric, sparse matrices with partial pivoting. It is also
applicable to rectangular systems of equations.

``umfpack``
-----------

The ``umfpack`` module interfaces the `UMFPACK
<http://www.cise.ufl.edu/research/sparse/umfpack>`_ factorization
package. UMFPACK computes the LU factorization of a general matrix with partial
pivoting. It is also applicable to rectangular systems.

The main difference between the ``superlu`` and ``umfpack`` modules resides in
the way the factorization is performed internally. SuperLU works with the
concept of *supernodes* leading naturally to parallelism in the
factorization. If available, a custom-built SuperLU library for multi-core
processors can be supplied to Pysparse in place of the default library. Both
factorization packages rely intensively on the BLAS to operate on dense
sub-matrices. Provided the BLAS library supplied to Pysparse was compiled with
multi-threading, some level of parallelism will also be available in
UMFPACK. However, a rough empirical observation is that UMFPACK is often faster
than SuperLU on mono-processor machines.

``jdsym``
---------

The ``jdsym`` module provides an implementation of the JDSYM algorithm, that is
conveniently callable from Python. The JDSYM algorithm computes solutions of
large sparse symmetric (genralised or standard) eigenvalue problems. JDSYM is an
implementation of the Jacobi-Davidson method, optimized for symmetric matrices.


Prerequisites
=============

- `NumPy <http://numpy.scipy.org>`_
- Optionally, a custom-built UMFPACK and/or SuperLU

Installing Pysparse
===================

``python setup.py install``
