.. Description of the jdsym module
.. module:: jdsym
.. _jdsym-page:

=================
Eigenvalue Solver
=================

The :mod:`jdsym` Module
=======================

The ``jdsym`` module provides an implementation of the JDSYM algorithm, that is
conveniently callable from Python. JDSYM is an eigenvalue solver to compute
eigenpairs of a generalised matrix eigenvalue problem of the form

.. math:: \mathbf{A} \mathbf{x} = \lambda \mathbf{M} \mathbf{x}
   :label: eq:python:2

or a standard eigenvalue problem of the form

.. math:: \mathbf{A} \mathbf{x} = \lambda \mathbf{x}
   :label: eq:python:3
  
where :math:`\mathbf{A}` is symmetric and :math:`\mathbf{M}` is symmetric
positive definite.

The module exports a single function:

.. function:: jdsym(A, M, K, kmax, tau, jdtol, itmax, linsolver, **kwargs)

   Implements Jacobi-Davidson iterative method to identify a given number of
   eigenvalues near a target value.

   :parameters:

        :A: the matrix :math:`\mathbf{A}` in :eq:`eq:python:2` or
             :eq:`eq:python:3`. ``A`` must provide the ``shape`` attribute and
             the ``matvec`` and ``matvec_transp`` methods.
        :M: the matrix :math:`\mathbf{M}`
             in :eq:`eq:python:2`. :math:`\mathbf{M}` must provide the ``shape``
             attribute and the ``matvec`` and ``matvec_transp`` methods. If the
             standard eigenvalue problem :eq:`eq:python:3` is to be solved,
             ``M`` should be set to ``None``.
        :K: a preconditioner object that supplies the ``shape`` attribute and
             the ``precon`` method. If no preconditioner is to used, then the
             ``None`` value can be passed for this parameter.
        :kmax: the number of eigenpairs to be computed.
        :tau: the target value :math:`\tau`. Eigenvalues in the vicinity
             of :math:`\tau` will be computed.
        :jdtol: the convergence tolerance for
             eigenpairs :math:`(\lambda,\mathbf{x})`. The converged eigenpairs
             have a residual :math:`\|\mathbf{A} \mathbf{x} - \lambda \mathbf{M}
             \mathbf{x}\|_2` less than ``jdtol``.
        :itmax: an integer that specifies the maximum number of Jacobi-Davidson
             iterations to perform.
        :linsolver: a function that implements an iterative method for solving
             linear systems of equations. The function ``linsolver`` is required
             to conform to the standards mentioned in :ref:`itsolvers-page`.

   :keywords:

        :jmax: the maximum dimension of the search subspace (default: 25).
        :jmin: the dimension of the search subspace after a restart (default:
             10).
        :blksize: the block size used in the JDSYM algorithm (default: 1).
        :blkwise: is an integer that affects the convergence criterion if
             ``blksize`` is larger than 1 (default: 0).
        :V0: a NumPy array of rank one or two. It specifies the initial
             search subspace (default: a randomly generated initial search
             subspace).
        :optype: is an integer specifying the operator type used in the
             correction equation. If ``optype=1``, the non-symmetric version is
             used. If ``optype=2``, the  symmetric version is used (default: 2).
        :linitmax: the maximum number of steps taken in the inner iteration
             (iterative linear solver) (default: 200).
        :eps_tr: the tracking parameter (default: 1.0e-3).
        :toldecay: is a float value that influences the dynamic adaptation of
             the stopping criterion of the inner iteration (default: 1.5).
        :clvl: verbosity level. The higher the ``clvl`` parameter, the more
             output is sent to the standard output. ``clvl=0`` produces no
             output (default: 0).
        :strategy: is an integer specifying shifting and sorting strategy of
             JDSYM. ``strategy=0`` enables the default JDSYM
             algorithm. ``strategy=1`` enables JDSYM to avoid convergence to
             eigenvalues smaller than :math:`\tau` (default: 0).
        :projector: is used to keep the search subspace and the eigenvectors
             in a certain subspace. The parameter ``projector`` can be any
             Python object that has a ``shape`` attribute and a ``project``
             method. The ``project`` method takes a vector (a rank-1 NumPy
             array) as its sole argument and projects that vector in-place. This
             parameter can be used to implement the DIRPROJ and SAUG methods
             (default: no projection).
   :returns:
        :kconv: the number of converged eigenpairs.
        :lambda: a rank-1 NumPy array containing the converged eigenvalues.
        :Q: a rank-2 NumPy array containing the converged eigenvectors. The
           i-th eigenvector is accessed by ``Q[:,i]``.
        :it: an integer indicating the number of Jacobi-Davidson steps
           (outer iteration steps) performed.



Example: Maxwell Problem
------------------------

.. todo:: Update the timings below.

.. warning:: The timings below are Roman's old benchmarks. We should run them
   again.

The following code illustrates the use of the ``jdsym`` module.  Two
matrices :math:`\mathbf{A}` and :math:`\mathbf{M}` are read from files. A Jacobi
preconditioner from :math:`\mathbf{A} - \tau\mathbf{M}` is built. Then the JDSYM
eigensolver is called, calculating 5 eigenvalues near 25.0 and the associated
eigenvalues to an accuracy of :math:`10^{-10}`.  We set ``strategy=1``
to avoid convergence to the high-dimensional null space of
(:math:`\mathbf{A}`, :math:`\mathbf{M}`)::

    from pysparse import spmatrix, itsolvers, jdsym, precon

    A = spmatrix.ll_mat_from_mtx('edge6x3x5_A.mtx')
    M = spmatrix.ll_mat_from_mtx('edge6x3x5_B.mtx')
    tau = 25.0

    Atau = A.copy()
    Atau.shift(-tau, M)
    K = precon.jacobi(Atau)

    A = A.to_sss(); M = M.to_sss()
    k_conv, lmbd, Q, it  = jdsym.jdsym(A, M, K, 5, tau,
                                       1e-10, 150, itsolvers.qmrs,
                                       jmin=5, jmax=10, clvl=1, strategy=1)

This code takes 33.71 seconds to compute the five eigenpairs.  A native
C version, using the same computational kernels, takes 35.64 for the same
task. We expected the Python version to be slower due to the overhead generated
when calling the matrix-vector multiplication and the preconditioner, but
surprisingly the Python code was even a bit faster.
