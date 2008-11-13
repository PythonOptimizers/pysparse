.. Description of the itsolvers module
.. module:: itsolvers
.. _itsolvers-page:

=================
Iterative Solvers
=================


The :mod:`itsolvers` Module
===========================

The ``itsolvers`` module provides a set of iterative methods for solving linear
systems of equations.

The iterative methods are callable like ordinary Python functions. All these
functions expect the same parameter list, and all function return values also
follow a common standard.

Any user-defined iterative solvers should also follow these conventions, since
other PySparse modules rely on them (e.g. the ``jdsym`` module;
see :ref:`jdsym-page`).

Let's illustrate the calling conventions, using the PCG method.

.. function:: info, iter, relres = pcg(A, b, x, tol, maxit[, K])

   Solve a linear system ``A x = b`` with the preconditioned conjugate gradient
   algorithm.

   :parameters:

        :A: The coefficient matrix of the linear system of equations. ``A``
             must provide the ``shape`` attribute and the ``matvec`` and
             ``matvec_transp`` methods for multiplying with a vector.
        :b: The right-hand-side of the linear system as a rank-1 NumPy array.
        :x: A rank-1 NumPy array. Upon entry, ``x`` holds the initial guess. On
             exit, ``x`` holds an approximate solution of the linear system.
        :tol: A float value representing the requested error tolerance. The
             exact meaning of this parameter depends on the actual iterative
             solver.
        :maxit: An integer that specifies the maximum number of iterations to
             be executed.
        :K: A preconditioner object that supplies the ``shape`` attribute and
             the ``precon`` method.

   :returns:

        :info: an integer that contains the exit status of the iterative solver.
           ``info >= 0`` indicates, that ``x`` holds an acceptable solution, and
           ``info < 0`` indicates an error condition. ``info`` has one of the
           following values:

           +2. iteration converged, residual is as small as seems reasonable on this machine,

           +1. iteration converged, ``b = 0``, so the exact solution is ``x = 0``.

           +0. iteration converged, relative error appears to be less than ``tol``.

           -1. iteration did not converge, maximum number of iterations was reached.

           -2. iteration did not converge, the system involving the preconditioner was ill-conditioned.

           -3. iteration did not converge, an inner product of the form :math:`\mathbf{x}^T \mathbf{K}^{-1} \mathbf{x}` was not positive, so the preconditioning matrix :math:`\mathbf{K}` does not appear to be positive definite.

           -4. iteration did not converge, the matrix ``A`` appears to be very ill-conditioned.

           -5. iteration did not converge, the method stagnated.

           -6. iteration did not converge, a scalar quantity became too small or too large to continue computing.

        :iter: the number of iterations performed.

        :relres: the relative residual at the approximate solution computed by
           the iterative method. What this actually is depends on the actual
           iterative method used.

The iterative solvers may accept additional parameters, which are passed as
keyword arguments.

Note that not all iterative solvers check for all above error conditions.


``itsolvers`` Module Functions
------------------------------

The module functions defined in the ``precon`` module implement
various iterative methods (PCG, MINRES, QMRS and CGS). The parameters and return
values conform to the conventions described above.

.. function:: info, iter, relres = pcg(A, b, x, tol, maxit[, K])

   Implements the Preconditioned Conjugate Gradient method.

.. function:: info, iter, relres = minres(A, b, x, tol, maxit[, K])

   Implements the MINRES method.

.. function:: info, iter, relres = qmrs(A, b, x, tol, maxit[, K])

   Implements the QMRS method.

.. function:: info, iter, relres = cgs(A, b, x, tol, maxit[, K])

   Implements the CGS method.


Example: Solving the Poisson System
-----------------------------------

Let's solve the Poisson system

.. math:: \mathbf{L} \mathbf{x} = \mathbf{1},
   :label: eq:python:1
  
using the PCG method. :math:`\mathbf{L}` is the 2D Poisson matrix, introduced in
:ref:`spmatrix-page`, and :math:`\mathbf{1}` is a vector with all
entries equal to one.

The Python solution for this task looks as follows::

    from pysparse import spmatrix, precon, itsolvers
    import numpy
    n = 300
    L = poisson2d_sym_blk(n)
    b = numpy.ones(n*n)
    x = numpy.empty(n*n)
    info, iter, relres = itsolvers.pcg(L.to_sss(), b, x, 1e-12, 2000)

The code makes use of the Python function ``poisson2d_sym_blk``,
which was defined in :ref:`spmatrix-page`.

Incorporating e.g. a SSOR preconditioner is straightforward::

    from pysparse import spmatrix, precon, itsolvers
    import numpy
    n = 300
    L = poisson2d_sym_blk(n)
    b = numpy.ones(n*n)
    x = numpy.empty(n*n)
    S = L.to_sss()
    Kssor = precon.ssor(S)
    info, iter, relres = itsolvers.pcg(S, b, x, 1e-12, 2000, Kssor)

The Matlab solution (without preconditioner) may look as follows:

.. code-block:: matlab

   n = 300;
   L = poisson2d_kron(n);
   [x,flag,relres,iter] = pcg(L, ones(n*n,1), 1e-12, 2000, ...
                              [], [], zeros(n*n,1));


Performance comparison with Matlab and native C
-----------------------------------------------

To evaluate the performance of the Python implementation we solve the 2D Poisson
system :eq:`eq:python:1` using the PCG method. The Python timings are compared
with results of a Matlab and a native C implementation.

The native C and the Python implementation use the same core algorithms for PCG
method and the matrix-vector multiplication. On the other hand, C reads the
matrix from an external file instead of building it on the fly. In contrast to
the Python implementation, the native C version does not suffer from the
overhead generated by the runtime argument parsing and calling overhead.

.. _python-vs-matlab-vs-c:

   **Table.** Performance comparison of Python, Matlab and native
   C implementations to solve the linear system :eq:`eq:python:1` without
   preconditioning. The execution times are given in seconds. *Assembly* is the
   time for constructing the matrix (or reading it from a file in the case of
   native C).  *Solve* is the time spent in the PCG solver. *Total* is the sum
   of *Assembly* and *Solve*. Matlab version 6.0 Release 12 was used for these
   timings.

   +----------+-------+----------+---------+--------+
   | Function | Size  | Assembly | Solve   | Total  |
   +----------+-------+----------+---------+--------+
   | Python   | n=100 | 0.03     | 1.12    | 1.15   |
   +----------+-------+----------+---------+--------+
   |          | n=300 | 0.21     | 49.65   | 49.86  |
   +----------+-------+----------+---------+--------+
   |          | n=500 | 0.62     | 299.39  | 300.01 |
   +----------+-------+----------+---------+--------+
   | Native C | n=100 | 0.30     | 0.96    | 1.26   |
   +----------+-------+----------+---------+--------+
   |          | n=300 | 3.14     | 48.38   | 51.52  |
   +----------+-------+----------+---------+--------+
   |          | n=500 | 10.86    | 288.67  | 299.53 |
   +----------+-------+----------+---------+--------+
   | Matlab   | n=100 | 0.21     | 8.85    | 9.06   |
   +----------+-------+----------+---------+--------+
   |          | n=300 | 2.05     | 387.26  | 389.31 |
   +----------+-------+----------+---------+--------+
   |          | n=500 | 6.23     | 1905.67 | 1911.8 |
   +----------+-------+----------+---------+--------+

This table shows the execution times for the Python, the
Matlab and the native C implementation for solving the linear system
:eq:`eq:python:1`. Matlab is not only slower when building the matrix, also
the matrix-vector multiplication seems to be implemented
inefficiently. Considering *Solve*, the performance of Python and
native C is comparable. The Python overhead is under a factor of 4.
