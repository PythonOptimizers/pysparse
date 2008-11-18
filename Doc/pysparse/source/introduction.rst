.. Introduction to Pysparse

========================
Introduction to Pysparse
========================

PySparse extends the Python interpreter by a set of sparse matrix types holding
double precision values. PySparse also includes modules that implement

- Iterative Krylov methods for solving linear systems of equations,
- Diagonal (Jacobi) and SSOR preconditioners,
- Interfaces to direct solvers for sparse linear systems of equations (SuperLU
  and UMFPACK),
- A Jacobi-Davidson eigenvalue solver for the symmetric, generalised matrix
  eigenvalue problem (JDSYM),
- Low-level C classes to represent and manipulate sparse matrices,
- High-level Python classes with operator overloading to perform usual
  operations on matrices.

Most of the above modules are implemented as C extension modules for maximum
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


Testing Pysparse
================

From the ``Test`` directory, ``testSuperLU`` runs a series of tests to exercise
the various options of the SuperLU direct solver::

    $ python testSuperlu.py
	      Test    RelErr       Tol    nnz(A)  nnz(L+U)    Fact   Solve
	-------------------------------------------------------------------
	poi1d-dflt  1.78e-12  2.25e-07     99999    199998    0.06    0.00
    .	poi1d-size  1.78e-12  2.25e-07     99999    199998    0.05    0.00
    .	poi1d-relx  1.78e-12  2.25e-07     99999    200758    0.06    0.00
    .	poi1d-trsh  1.78e-12  2.25e-07     99999    199998    0.06    0.00
    .	poi1d-prm0  1.44e-12  2.25e-07     99999    199998    0.04    0.00
    .	poi1d-prm1  1.78e-12  2.25e-07     99999    199998    0.07    0.00
    .	poi1d-prm2  1.78e-12  2.25e-07     99999    199998    0.06    0.00
    .	poi1d-prm3  1.44e-12  2.25e-07     99999    200000    0.06    0.00
    .	poi2d-dftl  2.55e-16  3.60e-12    119600   1952434    0.47    0.02
    .	poi2d-size  3.39e-16  3.60e-12    119600   1952434    0.42    0.02
    .	poi2d-relx  2.60e-16  3.60e-12    119600   2000252    0.48    0.02
    .	poi2d-trsh  2.55e-16  3.60e-12    119600   1952434    0.48    0.02
    .	poi2d-prm0  2.69e-15  3.60e-12    119600  16000398    5.48    0.11
    .	poi2d-prm1  7.00e-16  3.60e-12    119600   3506336    0.89    0.03
    .	poi2d-prm2  2.55e-16  3.60e-12    119600   1952434    0.47    0.02
    .	poi2d-prm3  2.24e-15  3.60e-12    119600   3472176    0.82    0.03
    .	spdgs-trsh  4.44e-16  2.22e-14     29998     39998    0.01    0.00
    .	spdgs-prm0  4.44e-16  2.22e-14     29998     40000    0.01    0.00
    .	spdgs-prm1  4.44e-16  2.22e-14     29998     40002    0.01    0.00
    .	spdgs-prm3  4.44e-16  2.22e-14     29998     40002    0.01    0.00
    .
    ----------------------------------------------------------------------
    Ran 20 tests in 12.675s

    OK

There is a corresponding test script for UMFPACK, ``testUmfpack``::

    $python testUmfpack.py
	  RelErr       Tol    nnz(A)    nnz(L)    nnz(U)    Fact   Solve
	-----------------------------------------------------------------
	1.44e-12  2.25e-07    149998     99999     99999    0.13    0.01
    .	1.44e-12  2.25e-07    149998     99999     99999    0.12    0.01
    .	1.44e-12  2.25e-07    149998     99999     99999    0.12    0.01
    .	1.44e-12  2.25e-07    149998     99999     99999    0.12    0.01
    .	1.44e-12  2.25e-07    149998     99999     99999    0.12    0.01
    .	1.44e-12  2.25e-07    149998     99999     99999    0.14    0.01
    .	1.55e-17  3.60e-12    199200   1081911   1081911    0.54    0.03
    .	1.55e-17  3.60e-12    199200   1081911   1081911    0.53    0.03
    .	2.50e-17  3.60e-12    199200   1081911   1081911    0.53    0.03
    .	2.50e-17  3.60e-12    199200   1081911   1081911    0.53    0.03
    .	1.55e-17  3.60e-12    199200   1081911   1081911    0.53    0.03
    .	1.64e-17  3.60e-12    199200   1489438   2166768    1.00    0.04
    .	4.44e-16  2.22e-14     29998     19999     19999    0.03    0.00
    .	4.44e-16  2.22e-14     29998     19999     19999    0.02    0.00
    .	4.44e-16  2.22e-14     29998     19999     19999    0.02    0.00
    .	4.44e-16  2.22e-14     29998     19999     19999    0.02    0.00
    .	4.44e-16  2.22e-14     29998     19999     19999    0.03    0.00
    .	4.44e-16  2.22e-14     29998     19999     19999    0.02    0.00
    .
    ----------------------------------------------------------------------
    Ran 18 tests in 8.486s

