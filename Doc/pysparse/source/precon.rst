.. Description of the precon module
.. module:: precon

===============
Preconditioners
===============

The :mod:`precon` Module
========================

The ``precon`` module provides preconditioners, which can be used e.g.\ for the
iterative methods implemented in the in the ``itsolvers`` module or the JDSYM
eigensolver (in the ``jdsym`` module).

In the Pysparse framework, any Python object that has the following
properties can be used as a preconditioner:

- a ``shape`` attribute, which returns a 2-tuple describing the dimension of the
  preconditioner,
- a ``precon`` method, that accepts two vectors :math:`\mathbf{x}`
  and :math:`\mathbf{y}`, and applies the preconditioner to :math:`\mathbf{x}`
  and stores the result in :math:`\mathbf{y}`. Both :math:`\mathbf{x}`
  and :math:`\mathbf{y}` are double precision, rank-1 NumPy arrays of
  appropriate size.

The ``precon`` module implements two new object types ``jacobi`` and ``ssor``,
representing Jacobi and SSOR preconditioners.


``precon`` Module Functions
---------------------------

.. function:: jacobi(A, omega=1.0, steps=1)

   Creates a ``jacobi`` object, representing the Jacobi preconditioner. The
   parameter :math:`\mathbf{A}` is the system matrix used for the Jacobi
   iteration. The matrix needs to be subscriptable using two-dimensional
   indices, so e.g.\ an ``ll_mat`` object would work.  The optional
   parameter :math:`\omega`, which defaults to 1.0, is the weight parameter.
   The optional ``steps`` parameter (defaults to 1) specifies the number of
   iteration steps.

.. function:: ssor(A, omega=1.0, steps=1)

   Creates a ``ssor`` object, representing the SSOR preconditioner.  The
   parameter :math:`\mathbf{A}` is the system matrix used for the SSOR
   iteration. The matrix :math:`\mathbf{A}` has to be an object of type
   ``sss_mat``.  The optional parameter :math:`\omega`, which defaults to 1.0,
   is the relaxation parameter. The optional ``steps`` parameter (defaults to 1)
   specifies the number of iteration steps.


``jacobi`` and ``ssor`` Objects
-------------------------------

Both ``jacobi`` and ``ssor`` objects provide the ``shape`` attribute and the
``precon`` method, that every preconditioner object in the PySparse framework
must implement.  Apart from that, there is nothing noteworthy to say about these
objects.

Example: Diagonal Preconditioner
--------------------------------

The diagonal preconditioner is just a special case of the Jacobi preconditioner,
with ``omega=1.0`` and ``steps=1``, which happen to be the default values of
these parameters.

It is however easy to implement the diagonal preconditioner using a
Python class::

       class diag_prec:
             def __init__(self, A):
                 self.shape = A.shape
                 n = self.shape[0]
                 self.dinv = numpy.empty(n)
                 for i in xrange(n):
                     self.dinv[i] = 1.0 / A[i,i]
             def precon(self, x, y):
                 numpy.multiply(x, self.dinv, y)

So::

    >>> D1 = precon.jacobi(A, 1.0, 1)

and::

    >>> D2 = diag_prec(A)

yield functionally equivalent preconditioners. ``D1`` is probably faster than
``D2``, because it is fully implemented in C.
