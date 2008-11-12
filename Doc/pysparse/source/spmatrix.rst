.. Description of the C-level sparse matrix types in spmatrix
.. module:: spmatrix
.. _spmatrix-page:

=============================
Low-Level Sparse Matrix Types
=============================


The :mod:`spmatrix` Module
==========================

The ``spmatrix`` module is the foundation of the PySparse
package. It extends the Python interpreter by three new types named
``ll_mat``, ``csr_mat`` and ``sss_mat``. These types
represent sparse matrices in the LL-, the CSR- and SSS-formats
respectively (see :ref:`formats-page`). For all three
formats, double precision values (C type ``double``) are used to
represent the non-zero entries.

The common way to use the ``spmatrix`` module is to first build a
matrix in the LL-format. The LL-matrix is manipulated until it has its
final shape and content. Afterwards it may be converted to either the
CSR- or SSS-format, which needs less memory and allows for fast
matrix-vector multiplications.

A ``ll_mat`` object can be created from scratch, by reading data
from a file (in MatrixMarket format) or as a result of matrix
operation (as e.g.\ a matrix-matrix multiplication). The
``ll_mat`` object supports manipulating (reading, writing,
add-updating) single entries or sub-matrices.

``csr_mat`` and ``sss_mat`` are not constructed directly,
instead they are created by converting ``ll_mat`` objects. Once
created, ``csr_mat`` and ``sss_mat`` objects cannot be
manipulated. Their purpose is to support efficient matrix-vector
multiplications.


``spmatrix`` module functions
-----------------------------

.. function:: ll_mat(n, m, sizeHint=1000)

   Creates a ``ll_mat`` object, that represents a general, all
   zero :math:`m\times n` matrix. The optional ``sizeHint`` parameter specifies
   the number of non-zero entries for which space is allocated initially.

   If the total number of non-zero elements of the final matrix is known
   (approximately), this number can be passed as ``sizeHint``. This
   will avoid costly memory reallocations.

.. function:: ll_mat_sym(n, sizeHint=1000)

   Creates a ``ll_mat`` object, that represents a *symmetric*,
   all zero :math:`n \times n` matrix. The optional ``sizeHint`` parameter
   specifies, how much space is initially allocated for the matrix.

.. function:: ll_mat_from_mtx(fileName)

   Creates a ``ll_mat`` object from a file named ``fileName``,
   which must be in `MatrixMarket Coordinate format
   <http://math.nist.gov/MatrixMarket/formats.html>`_. Depending on the
   file content, either a symmetric or a general sparse matrix is
   generated.

.. function:: matrixmultiply(A, B)

   Computes the matrix-matrix multiplication :math:`\mathbf{C} := \mathbf{A}
   \mathbf{B}` and returns the result :math:`\mathbf{C}` as a new ``ll_mat``
   object representing a general sparse matrix. The
   parameters :math:`\mathbf{A}` and :math:`\mathbf{B}` are expected to be
   objects of type ``ll_mat``.

.. function:: dot(A, B)

   Computes the *dot-product* :math:`\mathbf{C} := \mathbf{A}^T \mathbf{B}` and
   returns the result :math:`\mathbf{C}` as a new ``ll_mat`` object representing
   a general sparse matrix. The parameters :math:`\mathbf{A}`
   and :math:`\mathbf{B}` are expected to be objects of type ``ll_mat``.


``ll_mat`` objects
------------------

``ll_mat`` objects represent matrices stored in the LL format, which are
described in :ref:`formats-page`. ``ll_mat`` objects come in two flavours:
general matrices and symmetric matrices.  For symmetric matrices only the
non-zero entries in the lower triangle are stored.  Write operations to the
strictly upper triangle are prohibited for the symmetric format.  The ``issym``
attribute of an ``ll_mat`` object can be queried to find out whether or not the
symmetric storage format is used.

The entries of a matrix can be accessed conveniently using two-dimensional array
indices. In the Python language, subscripts can be of any type (as it is
customary for dictionaries). A two-dimensional index can be regarded as
a 2-tuple (the brackets do not have to be written, so ``A[1,2]`` is the same as
``A[(1,2)]``). If both tuple elements are integers, then a single matrix element
is referenced.  If at least one of the tuple elements is a slice (which is also
a Python object), then a submatrix is referenced.

Subscripts have to be decoded at runtime. This task includes type checks,
extraction of indices from the 2-tuple, parsing of slice objects and index bound
checks.  Following Python conventions, indices start with 0 and wrap around
(so -1 is equivalent to the last index).

The following code creates an empty :math:`5 \times 5` matrix ``A``, sets
all diagonal elements to their respective row/column index and then
copies the value of ``A[0,0]`` to ``A[2,1]``::

    >>> from pysparse import spmatrix
    >>> A = spmatrix.ll_mat(5, 5)
    >>> for i in range(5):
    ...     A[i,i] = i+1
    >>> A[2,1] = A[0,0]
    >>> print A
    ll_mat(general, [5,5], [(0,0): 1, (1,1): 2, (2,1): 1,
    (2,2): 3, (3,3): 4, (4,4): 5])


The Python slice notation can be used to conveniently access sub-matrices.

    >>> print A[:2,:]     # the first two rows
    ll_mat(general, [2,5], [(0,0): 1, (1,1): 2])
    >>> print A[:,2:5]    # columns 2 to 4
    ll_mat(general, [5,3], [(2,0): 3, (3,1): 4, (4,2): 5])
    >>> print A[1:3,2:5]  # submatrix from row 1 col 2 to row 2 col 4
    ll_mat(general, [2,3], [(1,0): 3])

The slice operator always returns a new ``ll_mat`` object, containing a **copy**
of the selected submatrix.

Write operations to slices are also possible:

    >>> B = ll_mat(2, 2)           # create 2-by-2
    >>> B[0,0] = -1; B[1,1] = -1   # diagonal matrix
    >>> A[:2,:2] = B               # assign it to upper
    >>>                            # diagonal block of A
    >>> print A
    ll_mat(general, [5,5], [(0,0): -1, (1,1): -1, (2,1): 1, 
    (2,2): 3, (3,3): 4, (4,4): 5])

*TODO: Mention fancy indexing*


``ll_mat`` Object Attributes and Methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. class:: ll_mat

   .. attribute:: shape

      Returns a 2-tuple containing the shape of the matrix :math:`\mathbf{A}`,
      i.e. the number of rows and columns.

   .. attribute:: nnz

      Returns the number of non-zero entries stored in
      matrix :math:`\mathbf{A}`. If :math:`\mathbf{A}` is stored in symmetric
      format, only the number of non-zero entries in the lower triangle
      (including the diagonal) are returned.

   .. attribute:: issym

      Returns true (a non-zero integer) if matrix :math:`\mathbf{A}` is stored
      in the symmetric LL format, i.e. only the non-zero entries in the lower
      triangle are stored. Returns false (zero) if matrix :math:`\mathbf{A}` is
      stored in the general LL format.

   .. method:: matvec(x, y)

      Computes the sparse matrix-vector product :math:`\mathbf{y} := \mathbf{A}
      \mathbf{x}` where :math:`\mathbf{x}` and :math:`\mathbf{y}` are double
      precision, rank-1 NumPy arrays of appropriate size.

   .. method:: matvec_transp(x, y) 

      Computes the transposed sparse matrix-vector product :math:`\mathbf{y} :=
      \mathbf{A}^T \mathbf{x}` where :math:`\mathbf{x}` and :math:`\mathbf{y}`
      are double precision, rank-1 NumPy arrays of appropriate size. For
      ``sss_mat`` objects ``matvec_transp`` is equivalent to ``matvec``.

   .. method:: to_csr()

      This method converts a sparse matrix in linked list format to compressed
      sparse row format. Returns a newly allocated ``csr_mat`` object, which
      results from converting matrix :math:`\mathbf{A}`.

   .. method:: to_sss()

      This method converts a sparse matrix in linked list format to sparse
      skyline format.  Returns a newly allocated ``sss_mat`` object, which
      results from converting matrix :math:`\mathbf{A}`. This function works for
      ``ll_mat`` objects in both the symmetric and the general
      format. If :math:`\mathbf{A}` is stored in general format, only the
      entries in the lower triangle are used for the conversion. No check
      whether :math:`\mathbf{A}` is symmetric is performed.

   .. method:: export_mtx(fileName, precision=6)

      Exports the matrix :math:`\mathbf{A}` to file named ``fileName``. The
      matrix is stored in `MatrixMarket Coordinate format
      <http://math.nist.gov/MatrixMarket/formats.html>`_. Depending on the
      properties of the ``ll_mat`` object :math:`\mathbf{A}` the generated file
      either uses the symmetric or a general MatrixMarket Coordinate format.
      The optional parameter ``precision`` specifies the number of decimal
      digits that are used to express the non-zero entries in the output file.

   .. method:: shift(sigma, M)

      Performs the daxpy operation :math:`\mathbf{A} \leftarrow \mathbf{A} +
      \sigma \mathbf{M}`. The parameter :math:`\sigma` is expected to be
      a Python Float object. The parameter :math:`\mathbf{M}` is expected to an
      object of type ``ll_mat`` of compatible shape.

   .. method:: copy()

      Returns a new ``ll_mat`` object, that represents a copy of the ``ll_mat``
      object :math:`\mathbf{A}`. So::

          >>> B = A.copy()

      is equivalent to::

          >>> B = A[:,:]

      On the other hand::

          >>> B = A.copy()

      is *not* the same as::

          >>> B = A

      The latter version only returns a reference to the same object and assigns
      it to ``B``. Subsequent changes to ``A`` will therefore also be visible in
      ``B``.

   .. method:: update_add_mask(B, ind0, ind1, mask0, mask1)

      This method is provided for efficiently assembling global finite element
      matrices. The method adds the matrix :math:`\mathbf{B}` to entries of
      matrix :math:`\mathbf{A}`. The indices of the entries to be updated are
      specified by the integer arrays ``ind0`` and ``ind1``. The individual
      updates are enabled or disabled using the ``mask0`` and ``mask1`` arrays.

      The operation is equivalent to the following Python code::

          for i in range(len(ind0)):
              for j in range(len(ind1)):
                  if mask0[i] and mask1[j]:
                     A[ind0[i],ind1[j]] += B[i,j]

      All five parameters are NumPy arrays. :math:`\mathbf{B}` is a rank-2
      array. The four remaining parameters are rank-1 arrays. Their length
      corresponds to either the number of rows or the number of columns
      of :math:`\mathbf{B}`.

      This method is not supported for ``ll_mat`` objects of symmetric type,
      since it would generally result in an non-symmetric matrix.
      ``update_add_mask_sym`` must be used in that case.  Attempting to call
      this method using a ``ll_mat`` object of symmetric type will raise an
      exception.

   .. method:: update_add_mask_sym(B, ind, mask)

      This method is provided for efficiently assembling symmetric global finite
      element matrices. The method adds the matrix :math:`\mathbf{B}` to entries
      of matrix :math:`\mathbf{A}`. The indices of the entries to be updated are
      specified by the integer array ``ind``. The individual updates are enabled
      or disabled using the ``mask`` array.

      The operation is equivalent to the following Python code::

          for i in range(len(ind)):
              for j in range(len(ind)):
                  if mask[i]:
                     A[ind[i],ind[j]] += B[i,j]

      The three parameters are all NumPy arrays. :math:`\mathbf{B}` is a rank-2
      array representing a square matrix. The two remaining parameters are
      rank-1 arrays. Their length corresponds to the order of
      matrix :math:`\mathbf{B}`.

   .. method:: update_add_at(val, irow, jcol):

      Add in place the elements of the vector ``val`` at the indices given by
      the two arrays ``irow`` and ``jcol``. The operation is equivalent to::

          for i in range(len(b)):
              A[irow[i],jcol[i]] += val[i]

   .. method:: generalize()

      Convert ``ll_mat`` object to non-symmetric form in place.

   .. method:: compress()

      Frees memory by reclaiming unused space in the internal data structure.
      Returns the number of elements freed.

   .. method:: norm(p)

     Returns the ``p``-norm of a matrix, where ``p`` is a string. If ``p='1'``,
     the 1-norm is returned, if ``p='inf'``, the infinity-norm is returned, and
     if ``p='fro'``, the Frobenius norm is returned.

     .. note::
        The 1 and infinity norm are not implemented for symmetric matrices.

   .. method:: keys()

      Return a list of tuples (i,j) of the nonzero matrix entries.

   .. method:: values()

      Return a list of the nonzero matrix entries as floats.

   .. method:: items()

      Return a list of tuples (indices, value) of the nonzero entries keys and
      values. The indices are themselves tuples (i,j) of row and column values.

   .. method:: scale(sigma)

      Scale each element in the matrix by the constant sigma.

   .. method:: take(val, irow, jcol)

      Extract elements at positions ``(irow[i], jcol[i])`` and place them in the
      array ``val``. In other words::

           for i in range(len(val)): val[i] = A[irow[i],jcol[i]]

   .. method:: put(val, irow, jcol)

      ``put`` takes the opposite tack to ``take``. Place the values in ``val``
      at positions given by ``irow`` and ``jcol``::

           for i in range(len(val)): A[irow[i],jcol[i]] = val[i]

   .. method:: delete_rows(mask)

      Delete rows in place. If ``mask[i] == 0``, the ``i``-th row is
      deleted. This operation does not simple zero out rows, they are *removed*,
      i.e., the resulting matrix is *smaller*.

   .. method:: delete_cols(mask)

      Similar to ``delete_rows`` only with columns.

   .. method:: delete_rowcols(mask)

      If ``mask[i] == 0`` both the ``i``-th row and the ``i``-th column are
      deleted.

   .. method:: find()

      Returns a triple ``(val,irow,jcol)`` of Numpy arrays containing the matrix
      in coordinate format, there is a nonzero element with value ``val[i]`` in
      position ``(irow[i],jcol[i])``.


``csr_mat`` and ``sss_mat`` Objects
-----------------------------------

``csr_mat`` objects represent matrices stored in the CSR format, which are
described in :ref:`formats-page`.  ``sss_mat`` objects represent matrices stored
in the SSS format (c.f. :ref:`formats-page`).  The only way to create
a ``csr_mat`` or a ``sss_mat`` object is by conversion of a ``ll_mat`` object
using the ``to_csr()`` or the ``to_sss()`` method respectively. The purpose of
the ``csr_mat`` and the ``to_sss()`` objects is to provide fast matrix-vector
multiplications for sparse matrices. In addition, a matrix stored in the CSR or
SSS format uses less memory than the same matrix stored in the LL format, since
the ``link`` array is not needed.

``csr_mat`` and ``sss_mat`` objects do not support two-dimensional indices to
access matrix entries or sub-matrices.  Again, their purpose is to provide fast
matrix-vector multiplication.

``csr_mat`` and ``sss_mat`` Object Attributes and Methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. class:: csr_mat
.. class:: sss_mat

   .. attribute:: shape

      Returns a 2-tuple containing the shape of the matrix :math:`\mathbf{A}`,
      i.e. the number of rows and columns.

   .. attribute:: nnz

      Returns the number of non-zero entries stored in
      matrix :math:`\mathbf{A}`. If :math:`\mathbf{A}` is an ``sss_mat`` object,
      the non-zero entries in the strictly upper triangle are not counted.

   .. method:: matvec(x, y)

      Computes the sparse matrix-vector product :math:`\mathbf{y} := \mathbf{A}
      \mathbf{x}` where :math:`\mathbf{x}` and :math:`\mathbf{y}` are double
      precision, rank-1 NumPy arrays of appropriate size.

   .. method:: matvec_transp(x, y) 

      Computes the transposed sparse matrix-vector product :math:`\mathbf{y} :=
      \mathbf{A}^T \mathbf{x}` where :math:`\mathbf{x}` and :math:`\mathbf{y}`
      are double precision, rank-1 NumPy arrays of appropriate size. For
      ``sss_mat`` objects ``matvec_transp`` is equivalent to ``matvec``.


Example: 2D-Poisson matrix
==========================

This section illustrates the use of the ``spmatrix`` module to
build the well known 2D-Poisson matrix resulting from a :math:`n \times n`
square grid.::

       def poisson2d(n):
           n2 = n*n
           L = spmatrix.ll_mat(n2, n2, 5*n2-4*n)
           for i in range(n):
               for j in range(n):
                   k = i + n*j
                   L[k,k] = 4
                   if i > 0:
                      L[k,k-1] = -1
                   if i < n-1:
                      L[k,k+1] = -1
                   if j > 0:
                      L[k,k-n] = -1
                   if j < n-1:
                      L[k,k+n] = -1
       return L

Using the symmetric variant of the ``ll_mat`` object, this gets
even shorter.::

     def poisson2d_sym(n):
         n2 = n*n
         L = spmatrix.ll_mat_sym(n2, 3*n2-2*n)
         for i in range(n):
             for j in range(n):
                 k = i + n*j
                 L[k,k] = 4
                 if i > 0:
                    L[k,k-1] = -1
                 if j > 0:
                    L[k,k-n] = -1
     return L

To illustrate the use of the slice notation to address sub-matrices,
let's build the 2D Poisson matrix using the diagonal and off-diagonal
blocks.::

        def poisson2d_sym_blk(n):
            n2 = n*n
            L = spmatrix.ll_mat_sym(n2, 2*n2-2*n)
            I = spmatrix.ll_mat_sym(n, n)
            P = spmatrix.ll_mat_sym(n, 2*n-1)
            for i in range(n):
                I[i,i] = -1
            for i in range(n):
                P[i,i] = 4
            if i > 0: P[i,i-1] = -1
            for i in range(0, n*n, n):
                L[i:i+n,i:i+n] = P
                if i > 0: L[i:i+n,i-n:i] = I
        return L


Performance comparison with Matlab
==================================

.. note ::
   These are Roman Gueus' tests. I am not sure on which type of machine they
   were performed. In the next section, we implement the above Poisson
   constructors by taking advantage of vectorization and run updated
   comparisons.

Let's compare the performance of three python codes above with the
following Matlab functions:

The Matlab function ``poisson2d`` is equivalent to the Python
function with the same name

.. code-block:: matlab

   function L = poisson2d(n)
            L = sparse(n*n);
            for i = 1:n
                for j = 1:n
                    k = i + n*(j-1);
                    L(k,k) = 4;
                    if i > 1, L(k,k-1) = -1; end
                    if i < n, L(k,k+1) = -1; end
                    if j > 1, L(k,k-n) = -1; end
                    if j < n, L(k,k+n) = -1; end
                end
            end

The function ``poisson2d_blk`` is an adaption of the Python function
``poisson2d_sym_blk`` (except for exploiting the symmetry, which is not
directly supported in Matlab).

.. code-block:: matlab

   function L = poisson2d_blk(n)
            e = ones(n,1);
            P = spdiags([-e 4*e -e], [-1 0 1], n, n);
            I = -speye(n);
            L = sparse(n*n);
            for i = 1:n:n*n
                L(i:i+n-1,i:i+n-1) = P;
                if i > 1, L(i:i+n-1,i-n:i-1) = I; end
                if i < n*n - n, L(i:i+n-1,i+n:i+2*n-1) = I; end
            end

The function ``poisson2d_kron`` demonstrates one of the most efficient 
ways to generate the 2D Poisson matrix in Matlab.

.. code-block:: matlab

   function L = poisson2d_kron(n)
            e = ones(n,1);
            P = spdiags([-e 2*e -e], [-1 0 1], n, n);
            L = kron(P, speye(n)) + kron(speye(n), P);

.. _python-vs-matlab:

   **Table.** Performance comparison of Python and Matlab functions to generate
   the 2D Poisson matrix. The execution times are given in seconds. Matlab
   version 6.0 Release 12 was used for these timings.

+------------------------------+-------+--------+---------+----------------+
| Function                     | n=100 |  n=300 |   n=500 | n=1000         |
+------------------------------+-------+--------+---------+----------------+
| Python ``poisson2d``         |  0.44 |   4.11 |   11.34 |  45.50         |
+------------------------------+-------+--------+---------+----------------+
| Python ``poisson2d_sym``     |  0.26 |   2.34 |   6.55  |  26.33         |
+------------------------------+-------+--------+---------+----------------+
| Python ``poisson2d_sym_blk`` |  0.03 |   0.21 |   0.62  |  2.22          |
+------------------------------+-------+--------+---------+----------------+
| Matlab ``poisson2d``         | 28.19 | 3464.9 | 38859.0 | :math:`\infty` |
+------------------------------+-------+--------+---------+----------------+
| Matlab ``poisson2d_blk``     | 6.85  | 309.20 | 1912.1  | :math:`\infty` |
+------------------------------+-------+--------+---------+----------------+
| Matlab ``poisson2d_kron``    | 0.21  | 2.05   | 6.23    | 29.96          |
+------------------------------+-------+--------+---------+----------------+

The execution times reported in Table python-vs-matlab_ clearly
show, that the Python implementation is superior to the Matlab
implementation. If the fastest versions are compared for both languages, Python
is approximately 10 times faster. Comparing the straight forward ``poisson2d``
versions, one is struck by the result that, the Matlab function is incredibly
slow. The Python version is more than three orders of magnitude faster!

.. This result really raises the doubt, whether Matlab's sparse matrix format is
.. appropriately chosen.

The performance difference between Python's ``poisson2d_sym`` and
  ``poisson2d_sym_blk`` indicates, that a lot of time is spent parsing indices.


Vectorization
=============

The ``put`` method of ``ll_mat`` objects allows us to operate on entire arrays
at a time. This is advantageous because the loop over the elements of an array
is performed at C level instead of in the Python script. For illustration, let's
rewrite the ``poisson2d``, ``poisson2d_sym`` and ``poisson2d_sym_blk``
constructors.

The ``put`` method can be used in ``poisson1d`` as so::

    from pysparse import spmatrix
    import numpy

    def poisson2d_vec(n):
        n2 = n*n
        L = spmatrix.ll_mat(n2, n2, 5*n2-4*n)
        e = numpy.ones(n)
        d = numpy.arange(n, dtype=numpy.int)
        din = d
        for i in xrange(n):
            # Diagonal blocks
            L.put(4*e, din, din)
            L.put(-e[1:], din[1:], din[:-1])
            L.put(-e[1:], din[:-1], din[1:])
            # Outer blocks
            L.put(-e, n+din, din)
            L.put(-e, din, n+din)
            din = d + i*n
        # Last diagonal block
        L.put(4*e, din, din)
        L.put(-e[1:], din[1:], din[:-1])
        L.put(-e[1:], din[:-1], din[1:])
        return L

And similarly in the symmetric version::

    def poisson2d_sym_vec(n):
        n2 = n*n
        L = spmatrix.ll_mat_sym(n2, 3*n2-2*n)
        e = numpy.ones(n)
        d = numpy.arange(n, dtype=numpy.int)
        din = d
        for i in xrange(n):
            # Diagonal blocks
            L.put(4*e, din, din)
            L.put(-e[1:], din[1:], din[:-1])
            # Outer blocks
            L.put(-e, n+din, din)
            din = d + i*n
        # Last diagonal block
        L.put(4*e, din, din)
        L.put(-e[1:], din[1:], din[:-1])
        return L

The time differences to construct matrices with and without vectorization can be
dramatic:

.. sourcecode:: ipython

   In [1]: from pysparse import poisson
   In [2]: import poisson_vec

   In [3]: %timeit -n10 -r3 L = poisson.poisson2d(100)
   10 loops, best of 3: 38.2 ms per loop
   In [4]: %timeit -n10 -r3 L = poisson_vec.poisson2d_vec(100)
   10 loops, best of 3: 8.11 ms per loop

   In [5]: %timeit -n10 -r3 L = poisson.poisson2d(300)
   10 loops, best of 3: 352 ms per loop
   In [6]: %timeit -n10 -r3 L = poisson_vec.poisson2d_vec(300)
   10 loops, best of 3: 44.8 ms per loop

   In [7]: %timeit -n10 -r3 L = poisson.poisson2d(500)
   10 loops, best of 3: 980 ms per loop
   In [8]: %timeit -n10 -r3 L = poisson_vec.poisson2d_vec(500)
   10 loops, best of 3: 110 ms per loop

   In [9]: %timeit -n10 -r3 L = poisson.poisson2d(1000)
   10 loops, best of 3: 4.02 s per loop
   In [10]: %timeit -n10 -r3 L = poisson_vec.poisson2d_vec(1000)
   10 loops, best of 3: 398 ms per loop
   
and for the symmetric versions:

.. sourcecode:: ipython

   In [18]: %timeit -n10 -r3 L = poisson.poisson2d_sym(100)
   10 loops, best of 3: 22.6 ms per loop
   In [19]: %timeit -n10 -r3 L = poisson_vec.poisson2d_sym_vec(100)
   10 loops, best of 3: 5.05 ms per loop

   In [20]: %timeit -n10 -r3 L = poisson.poisson2d_sym(300)
   10 loops, best of 3: 202 ms per loop
   In [21]: %timeit -n10 -r3 L = poisson_vec.poisson2d_sym_vec(300)
   10 loops, best of 3: 27 ms per loop

   In [22]: %timeit -n10 -r3 L = poisson.poisson2d_sym(500)
   10 loops, best of 3: 561 ms per loop
   In [23]: %timeit -n10 -r3 L = poisson_vec.poisson2d_sym_vec(500)
   10 loops, best of 3: 63.7 ms per loop

   In [24]: %timeit -n10 -r3 L = poisson.poisson2d_sym(1000)
   10 loops, best of 3: 2.26 s per loop
   In [25]: %timeit -n10 -r3 L = poisson_vec.poisson2d_sym_vec(1000)
   10 loops, best of 3: 224 ms per loop

From these numbers, it is obvious that vectorizing is crucial, especially for
large matrices. The gain in terms of time seems to be a factor of at least four
or five. Note that the last system has order one million.

Finally, the block version could be written as::

    def poisson2d_sym_blk_vec(n):
        n2 = n*n
        L = spmatrix.ll_mat_sym(n2, 3*n2-2*n)
        D = spmatrix.ll_mat_sym(n, 2*n-1)
        e = numpy.ones(n)
        d = numpy.arange(n, dtype=numpy.int)
        D.put(4*e, d, d)
        D.put(-e[1:], d[1:], d[:-1])
        P = spmatrix.ll_mat(n, n, n-1)
        P.put(-e,d,d)
        for i in xrange(n-1):
            L[i*n:(i+1)*n, i*n:(i+1)*n] = D
            L[(i+1)*n:(i+2)*n, i*n:(i+1)*n] = P
        # Last diagonal block
        L[n2-n:n2, n2-n:n2] = D
        return L

Here, ``put`` is sufficiently efficient that the benefit of constructing the
matrix by blocks is not apparent anymore. The slicing and block notation can
nevertheless be used for clarity. It could also be implemented as a combination
of ``find`` and ``put``, at the expense of memory consumption.

.. sourcecode:: ipython

   In [9]: %timeit -n10 -r3 L = poisson.poisson2d_sym_blk(1000)
   10 loops, best of 3: 246 ms per loop
   In [10]: %timeit -n10 -r3 L = poisson_vec.poisson2d_sym_blk_vec(1000)
   10 loops, best of 3: 246 ms per loop

The two best timings coinciding is a pure coincidence.
