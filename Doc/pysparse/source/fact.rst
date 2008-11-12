.. Describe the direct solvers
.. module:: superlu
.. _fact-page:

==============
Direct Solvers
==============

The Low-Level C Modules
=======================

The :mod:`superlu` Module
-------------------------

The ``superlu`` module interfaces the SuperLU library to make it usable by
Python code. SuperLU is a software package written in C, that is able to compute
an LU-factorisation of a general non-symmetric sparse matrix with
partial pivoting.

The ``superlu`` module exports a single function, called ``factorize``.

.. function:: factorize(A, **kwargs)

   The factorize function computes an LU-factorisation of the matrix ``A``.

   :parameters:

        :A: A ``csr_mat`` object that represents the matrix to be factorized.

   :keywords:

        :diag_pivot_thresh: the partial pivoting threshold, in the
                            interval :math:`[0,1]`. ``diag_pivot_thresh=0``
                            corresponds to no pivoting. ``diag_pivot_thresh=1``
                            corresponds to partial pivoting (default: 1.0).
        :drop_tol: the drop tolerance, in the
                            interval :math:`[0,1]`. ``drop_tol=0`` corresponds
                            to the exact factorization (default: 0.0).
        :relax: the degree of relaxing supernodes (default: 1).
        :panel_size: the maximum number of columns that form a panel (default:
                            10).
        :permc_spec: the matrix ordering used to control sparsity of the
                     factors:

                     0. natural ordering
                     1. MMD applied to the structure
                        of :math:`\mathbf{A}^T\mathbf{A}`
                     2. MMD applied to the structure of
                        :math:`\mathbf{A}^T + \mathbf{A}`
                     3. COLAMD, approximate minimum degree column ordering

                     (default: 2).

   :rtype: an object of type :class:`superlu_context`. This object encapsulates
           the L and U factors of ``A`` (see below).

.. note::

   The ``drop_tol`` has no effect in SuperLU version 2.0 and below. In SuperLU
   version 3.0 and above, the default value of ``permc_spec`` is 3.


``superlu_context`` Object Attributes and Methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. class:: superlu_context

   An abstract encapsulation of the LU factorization of a matrix by SuperLU.

   .. attribute:: shape 

      A 2-tuple describing the dimension of the  matrix factorized. It is equal
      to ``A.shape``.

   .. attribute:: nnz

      The ``nnz`` attribute holds the total number of nonzero entries stored in
      both the L and U factors.


   .. method:: solve(b, x, trans)

      The ``solve`` method accepts two rank-1 NumPy arrays ``b`` and ``x`` of
      appropriate size and assigns the solution of the linear system
      :math:`\mathbf{A}\mathbf{x} = \mathbf{b}` to ``x``. If the optional
      parameter ``trans`` is set to the string ``'T'``, the transposed system
      :math:`\mathbf{A}^T\mathbf{x} = \mathbf{b}` is solved instead.


Example: 2D Poisson Matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^

Let's now solve the 2D Poisson system :math:`\mathbf{A} \mathbf{x} = \mathbf{1}`
using an LU factorization. Here, :math:`\mathbf{A}` is the 2D Poisson matrix,
introduced in :ref:`spmatrix-page` and :math:`\mathbf{1}` is a vector with
all entries equal to one.

The Python solution for this task looks as follows::

    from pysparse import spmatrix, superlu
    import numpy 
    n = 100
    A = poisson2d_sym_blk(n)
    b = numpy.ones(n*n)
    x = numpy.empty(n*n)
    LU = superlu.factorize(A.to_csr(), diag_pivot_thresh=0.0)
    LU.solve(b, x)

The code makes use of the Python function :func:`poisson2d_sym_blk`, which
was defined in :ref:`spmatrix-page`.


Example: An Incomplete LU Factorization Preconditioner
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. warning::
   SuperLU 3.0 and above accept a ``drop_tol`` argument although the source
   files mention that incomplete factorization is *not implemented*. Therefore,
   changing ``drop_tol`` has no effect on the factorization at the moment and we
   must wait for it to be implemented. In the meantime, we can still demonstrate
   in this section how to implement an incomplete factorization preconditioner
   in Pysparse, even though in the present situation, it will be a *complete*
   factorization preconditioner!

Versions of SuperLU above 3.0 accept the ``drop_tol`` argument that allows the
computation of incomplete factors, realizing a tradeoff between computational
cost and factor density. The following example show how to use an incomplete LU
factorization as a preconditioner in any of the iterative methods of the
``itsolvers`` module::

    from pysparse import poisson, superlu, itsolvers
    import numpy

    class ILU_Precon:
        """
        A preconditioner based on an
        incomplete LU factorization.

        Input: A matrix in CSR format.
        Keyword argument: Drop tolerance.
        """
        def __init__(self, A, drop=1.0e-3):
            self.LU = superlu.factorize(A, drop_tol=drop)
            self.shape = self.LU.shape

        def precon(self, x, y):
            self.LU.solve(x,y)


    n = 300
    A = poisson.poisson2d_sym_blk(n).to_csr()   # Convert right away
    b = numpy.ones(n*n)
    x = numpy.empty(n*n)

    K = ILU_Precon(A)
    info, niter, relres = itsolvers.pcg(A, b, x, 1e-12, 2000, K)

.. note::

   Note that the 2D Poisson matrix is symmetric and positive definite, although
   barely. Indeed its smallest eigenvalue is
   :math:`2 (1 - \cos(\pi/(n+1))) \approx (\pi/(n+1))^2`. Therefore, a Cholesky
   factorization would be more appropriate. In the future, we intend to
   interface the `Cholmod <http://www.cise.ufl.edu/research/sparse/cholmod/>`_
   library.


The :mod:`umfpack` Module
-------------------------

*TODO*

Higher-Level Python Interfaces
==============================

The Abstract :mod:`directSolver` Module
---------------------------------------

.. automodule:: directSolver

.. autoclass:: PysparseDirectSolver
   :show-inheritance: 
   :members: 
   :inherited-members: 
   :undoc-members:

The :mod:`pysparseSuperLU` Module: A Higher-Level SuperLU Interface
-------------------------------------------------------------------

.. automodule:: pysparseSuperLU

.. autoclass:: PysparseSuperLUSolver
   :show-inheritance: 
   :members: 
   :inherited-members: 
   :undoc-members:

Example: The 2D Poisson System with SuperLU
-------------------------------------------

The solution of a 2D Poisson system with :class:`PysparseSuperLUSolver` may look
like this::

     from pysparse.pysparseMatrix import PysparseMatrix
     from pysparse.pysparseSuperLU import PysparseSuperLUSolver

     from poisson_vec import poisson2d_sym_blk_vec
     from numpy import ones
     from numpy.linalg import norm

     n = 200
     A = PysparseMatrix( matrix=poisson2d_sym_blk_vec(n) )
     x_exact = ones(n*n)/n
     b = A * x_exact
     LU = PysparseSuperLUSolver(A)
     LU.solve(b)
     print 'Factorization time: ', LU.factorizationTime
     print 'Solution time: ', LU.solutionTime
     print 'Error: ', norm(LU.sol - x_exact)/norm(x_exact)

The above script produces the output::

    Factorization time:  0.494116
    Solution time:  0.017096
    Error: 2.099685128150953e-14

Note that this example uses the vectorized Poisson constructors
of :ref:`spmatrix-page`.


The :mod:`pysparseUmfpack` Module: A Higher-Level UMFPACK Interface
-------------------------------------------------------------------

.. automodule:: pysparseUmfpack

.. autoclass:: PysparseUmfpackSolver
   :show-inheritance: 
   :members: 
   :inherited-members: 
   :undoc-members:

Example: The 2D Poisson System with UMFPACK
-------------------------------------------

The solution of a 2D Poisson system with :class:`PysparseUmfpackSolver` may look
like this::

     from poisson_vec import poisson2d_sym_blk_vec
     from numpy import ones
     from numpy.linalg import norm

     n = 200
     A = PysparseMatrix( matrix=poisson2d_sym_blk_vec(n) )
     x_exact = ones(n*n)/n
     b = A * x_exact
     LU = PysparseUmfpackSolver(A)
     LU.solve(b)
     print 'Factorization time: ', LU.factorizationTime
     print 'Solution time: ', LU.solutionTime
     print 'Error: ', norm(LU.sol - x_exact)/norm(x_exact)

This script produces the output::

     Factorization time:  0.520043
     Solution time:  0.031086
     Error: 1.10998989668e-15
