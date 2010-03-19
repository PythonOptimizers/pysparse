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


:mod:`spmatrix` module functions
--------------------------------

.. function:: ll_mat(n, m, sizeHint=1000)

   Creates a ``ll_mat`` object, that represents a general, all
   zero :math:`m \times n` matrix. The optional ``sizeHint`` parameter specifies
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

   Computes the *dot-product* :math:`\mathbf{C} := \mathbf{A^T B}` and
   returns the result :math:`\mathbf{C}` as a new ``ll_mat`` object representing
   a general sparse matrix. The parameters :math:`\mathbf{A}`
   and :math:`\mathbf{B}` are expected to be objects of type ``ll_mat``.


:class:`ll_mat` objects
-----------------------

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


Fancy Indexing
^^^^^^^^^^^^^^

There is flexibility in the way submatrices of ``ll_mat`` objects can be
accessed. In particular, rows and columns can be permuted arbitrarily and
submatrices need not be composed of consecutive rows or indices. Let's look at
an example. Below, the :func:`poisson1d` function assembles a Poisson matrix. We
come back to Poisson matrices later in this section. ::

    >>> from pysparse import poisson
    >>> n = 5
    >>> A = poisson.poisson1d(n)
    >>> print A   # Original matrix
    ll_mat(general, [5,5]):
     2.000000 -1.000000  --------  --------  -------- 
    -1.000000  2.000000 -1.000000  --------  -------- 
     -------- -1.000000  2.000000 -1.000000  -------- 
     --------  -------- -1.000000  2.000000 -1.000000 
     --------  --------  -------- -1.000000  2.000000 

    >>> print A[n-1:1:-1,1:n-1]  # Rows 2 through n-1 in reverse order,
    >>>                          # second through one before last col
    ll_mat(general, [3,3]):
     --------  -------- -1.000000
     -------- -1.000000  2.000000
    -1.000000  2.000000 -1.000000

    >>> print A[::-1,:]  # Reverse row order
    ll_mat(general, [5,5]):
     --------  --------  -------- -1.000000  2.000000 
     --------  -------- -1.000000  2.000000 -1.000000 
     -------- -1.000000  2.000000 -1.000000  -------- 
    -1.000000  2.000000 -1.000000  --------  -------- 
     2.000000 -1.000000  --------  --------  -------- 

    >>> print A[:,::-1]  # Reverse col order (same as above b/c A is symmetric)
    ll_mat(general, [5,5]):
     --------  --------  -------- -1.000000  2.000000 
     --------  -------- -1.000000  2.000000 -1.000000 
     -------- -1.000000  2.000000 -1.000000  -------- 
    -1.000000  2.000000 -1.000000  --------  -------- 
     2.000000 -1.000000  --------  --------  -------- 

    >>> print A[::-1,::-1]  # Reverse row and col order (same as original matrix)
    ll_mat(general, [5,5]):
     2.000000 -1.000000  --------  --------  -------- 
    -1.000000  2.000000 -1.000000  --------  -------- 
     -------- -1.000000  2.000000 -1.000000  -------- 
     --------  -------- -1.000000  2.000000 -1.000000 
     --------  --------  -------- -1.000000  2.000000

    >>> print A[1:3,3:]   # Rows 1 and 2, cols 3 and up
    ll_mat(general, [2,2]):
     --------  -------- 
    -1.000000  -------- 

    >>> print A[::2,::2]  # Every other row and col
    ll_mat(general, [3,3]):
     2.000000  --------  -------- 
     --------  2.000000  -------- 
     --------  --------  2.000000 

Keep in mind that as always with Python slices, the final index is never
included. Note also that slicing always returns a general matrix. Even though it
might be symmetric, both triangles are stored. Finally, slicing should be
applied to general matrices. If applied to symmetric matrices, only a partial
result is returned.

Fancy indexing can also be done with Python lists::

    >>> print A[ [1,4,2,0], ::2]
    ll_mat(general, [4,3]):
    -1.000000 -1.000000  -------- 
     --------  --------  2.000000 
     --------  2.000000  -------- 
     2.000000  --------  -------- 
    >>> p = [1,4,2,0]
    >>> q = [0,2,4]
    >>> print A[p,q]
    ll_mat(general, [4,3]):
    -1.000000 -1.000000  -------- 
     --------  --------  2.000000 
     --------  2.000000  -------- 
     2.000000  --------  --------

or with integer Numpy arrays::

   >>> idx0 = numpy.array([1,4,2,0], dtype=numpy.int)
   >>> idx1 = numpy.array([0,2,4], dtype=numpy.int)
   >>> print A[idx0,idx1]
   ll_mat(general, [4,3]):
   -1.000000 -1.000000  -------- 
    --------  --------  2.000000 
    --------  2.000000  -------- 
    2.000000  --------  --------

Finally, fancy indexing can be used to assign the same numerical value to
a submatrix::

    >>> A[:3,:3] = 7   # Assign value 7.0 to a principal submatrix
    >>> print A
    ll_mat(general, [5,5]):
     7.000000  7.000000  7.000000  --------  -------- 
     7.000000  7.000000  7.000000  --------  -------- 
     7.000000  7.000000  7.000000 -1.000000  -------- 
     --------  -------- -1.000000  2.000000 -1.000000 
     --------  --------  -------- -1.000000  2.000000 

Notice however that although the slice ``[0:3]`` appears to amount to the
list ``[0,1,2]``, the assignments ``A[:3,:3]=7`` and
``A.put([7,7,7], [0,1,2], [0,1,2])`` produce **very different** results.

.. warning:: For large-scale matrices, fancy indexing is most efficient when
             both index sets have the same type: two Python slices or two Python
             lists. When the index sets have different types, index arrays are
             built internally and this results in a performance hit.


:class:`ll_mat` Object Attributes and Methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. class:: ll_mat

   A general sparse matrix class in linked-list format which also allows the
   representation of symmetric matrices. Only the lower triangle of a symmetric
   matrix is kept in memory for efficiency.

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

      Returns a new ``ll_mat`` object that is a (deep) copy of the ``ll_mat``
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

   .. method:: update_add_at(val, irow, jcol)

      Add in place the elements of the vector ``val`` at the indices given by
      the two arrays ``irow`` and ``jcol``. The operation is equivalent to::

          for i in range(len(val)):
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
        The 1 and infinity norm are not yet implemented for symmetric matrices.

   .. method:: keys()

      Return a list of tuples ``(i,j)`` of the nonzero matrix entries.

   .. method:: values()

      Return a list of the nonzero matrix entries as floats.

   .. method:: items()

      Return a list of tuples ``(indices, value)`` of the nonzero entries keys
      and values. The indices are themselves tuples ``(i,j)`` of row and column
      values.

   .. method:: scale(sigma)

      Scale each element in the matrix by the constant ``sigma``.

   .. method:: take(val, irow, jcol)

      Extract elements at positions ``(irow[i], jcol[i])`` and place them in the
      array ``val``. In other words::

           for i in range(len(val)): val[i] = A[irow[i],jcol[i]]

   .. method:: put(val, irow, jcol)

      ``put`` takes the opposite tack to ``take``. Place the values in ``val``
      at positions given by ``irow`` and ``jcol``::

           for i in range(len(val)): A[irow[i],jcol[i]] = val[i]

      Here, ``irow`` and ``jcol`` can be Python lists or integer Numpy arrays.
      If either ``irow`` or ``jcol`` is omitted, it is replaced with
      ``[0, 1, 2, ...]``. Similarly, ``val`` can be a Python list, an integer
      Numpy array or a single scalar. If ``val`` is a scalar ``v``, it has the
      same effect as if it were the constant list or array ``[v, v, ..., v]``.

   .. method:: delete_rows(mask)

      Delete rows in place. If ``mask[i] == 0``, the ``i``-th row is
      deleted. This operation does not simply zero out rows, they are *removed*,
      i.e., the resulting matrix is *smaller*.

   .. method:: delete_cols(mask)

      Similar to ``delete_rows`` only with columns.

   .. method:: delete_rowcols(mask)

      If ``mask[i] == 0`` both the ``i``-th row and the ``i``-th column are
      deleted.

   .. method:: find()

      Returns a triple ``(val,irow,jcol)`` of Numpy arrays containing the matrix
      in coordinate format. There is a nonzero element with value ``val[i]`` in
      position ``(irow[i],jcol[i])``.


:class:`csr_mat` and :class:`sss_mat` Objects
---------------------------------------------

``csr_mat`` objects represent matrices stored in the CSR format, which are
described in :ref:`formats-page`.  ``sss_mat`` objects represent matrices stored
in the SSS format (c.f. :ref:`formats-page`).  The only way to create
a ``csr_mat`` or a ``sss_mat`` object is by conversion of a ``ll_mat`` object
using the :meth:`to_csr` or the :meth:`to_sss` method respectively. The purpose
of the ``csr_mat`` and the ``sss_mat`` objects is to provide fast matrix-vector
multiplications for sparse matrices. In addition, a matrix stored in the CSR or
SSS format uses less memory than the same matrix stored in the LL format, since
the ``link`` array is not needed.

``csr_mat`` and ``sss_mat`` objects do not support two-dimensional indices to
access matrix entries or sub-matrices.  Again, their purpose is to provide fast
matrix-vector multiplication.

:class:`csr_mat` and :class:`sss_mat` Object Attributes and Methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. class:: csr_mat

   A general sparse matrix class in compressed sparse row format which also
   allows the representation of symmetric matrices. Only the lower triangle of
   a symmetric matrix is kept in memory for efficiency.


.. class:: sss_mat

   A general sparse matrix class in sparse skyline format which also allows the
   representation of symmetric matrices. Only the lower triangle of a symmetric
   matrix is kept in memory for efficiency.

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
square grid::

       from pysparse import spmatrix

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
even shorter::

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
blocks::

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


Vectorization
=============

The ``put`` method of ``ll_mat`` objects allows us to operate on entire arrays
at a time. This is advantageous because the loop over the elements of an array
is performed at C level instead of in the Python script. 

If you need to ``put`` the same value in many places, ``put`` lets you specify
this value as a floating-point number instead of an array, e.g.::

    A.put(4.0, range(n), range(n))

is perfectly equivalent to::

    A.put(4*numpy.ones(n), range(n), range(n))

Moreover, if the second index set is omitted, it defaults to ``range(n)`` where
``n`` is the appropriate matrix dimension. So the above is again perfectly
equivalent to::

    A.put(4.0, range(n))

For illustration, let's rewrite the ``poisson2d``, ``poisson2d_sym`` and
``poisson2d_sym_blk`` constructors.

The ``put`` method can be used in ``poisson2d`` as so::

    from pysparse import spmatrix
    import numpy

    def poisson2d_vec(n):
        n2 = n*n
        L = spmatrix.ll_mat(n2, n2, 5*n2-4*n)
        d = numpy.arange(n2, dtype=numpy.int)
        L.put(4.0, d)
        L.put(-1.0, d[:-n], d[n:])
        L.put(-1.0, d[n:], d[:-n])
        for i in xrange(n):
            di = d[i*n:(i+1)*n]
            L.put(-1.0, di[1:], di[:-1])
            L.put(-1.0, di[:-1], di[1:])
        return L


And similarly in the symmetric version::

    def poisson2d_sym_vec(n):
        n2 = n*n
        L = spmatrix.ll_mat_sym(n2, 3*n2-2*n)
        d = numpy.arange(n2, dtype=numpy.int)
        L.put(4.0, d)
        L.put(-1.0, d[n:], d[:-n])
        for i in xrange(n):
            di = d[i*n:(i+1)*n]
            L.put(-1.0, di[:-1], di[1:])
        return L

The time differences to construct matrices with and without vectorization can be
dramatic. The following timings were generated on a 2.4GHz Intel Core2 Duo
processor:

.. sourcecode:: ipython

   In [1]: from pysparse import poisson
   In [2]: import poisson_vec

   In [3]: %timeit -n10 -r3 L = poisson.poisson2d(100)
   10 loops, best of 3: 38.2 ms per loop
   In [4]: %timeit -n10 -r3 L = poisson_vec.poisson2d_vec(100)
   10 loops, best of 3: 4.26 ms per loop

   In [5]: %timeit -n10 -r3 L = poisson.poisson2d(300)
   10 loops, best of 3: 352 ms per loop
   In [6]: %timeit -n10 -r3 L = poisson_vec.poisson2d_vec(300)
   10 loops, best of 3: 31.7 ms per loop

   In [7]: %timeit -n10 -r3 L = poisson.poisson2d(500)
   10 loops, best of 3: 980 ms per loop
   In [8]: %timeit -n10 -r3 L = poisson_vec.poisson2d_vec(500)
   10 loops, best of 3: 86.4 ms per loop

   In [9]: %timeit -n10 -r3 L = poisson.poisson2d(1000)
   10 loops, best of 3: 4.02 s per loop
   In [10]: %timeit -n10 -r3 L = poisson_vec.poisson2d_vec(1000)
   10 loops, best of 3: 333 ms per loop
   
and for the symmetric versions:

.. sourcecode:: ipython

   In [18]: %timeit -n10 -r3 L = poisson.poisson2d_sym(100)
   10 loops, best of 3: 22.6 ms per loop
   In [19]: %timeit -n10 -r3 L = poisson_vec.poisson2d_sym_vec(100)
   10 loops, best of 3: 2.48 ms per loop

   In [20]: %timeit -n10 -r3 L = poisson.poisson2d_sym(300)
   10 loops, best of 3: 202 ms per loop
   In [21]: %timeit -n10 -r3 L = poisson_vec.poisson2d_sym_vec(300)
   10 loops, best of 3: 20 ms per loop

   In [22]: %timeit -n10 -r3 L = poisson.poisson2d_sym(500)
   10 loops, best of 3: 561 ms per loop
   In [23]: %timeit -n10 -r3 L = poisson_vec.poisson2d_sym_vec(500)
   10 loops, best of 3: 53.8 ms per loop

   In [24]: %timeit -n10 -r3 L = poisson.poisson2d_sym(1000)
   10 loops, best of 3: 2.26 s per loop
   In [25]: %timeit -n10 -r3 L = poisson_vec.poisson2d_sym_vec(1000)
   10 loops, best of 3: 205 ms per loop

From these numbers, it is obvious that vectorizing is crucial, especially for
large matrices. The gain in terms of time seems to be a factor of at least four
or five. Note that the last system has order one million.

Finally, the block version could be written as::

    def poisson2d_vec_sym_blk(n):
        n2 = n*n
        L = spmatrix.ll_mat_sym(n2, 3*n2-2*n)
        D = spmatrix.ll_mat_sym(n, 2*n-1)
        d = numpy.arange(n, dtype=numpy.int)
        D.put(4.0, d)
        D.put(-1.0, d[1:], d[:-1])
        P = spmatrix.ll_mat_sym(n, n-1)
        P.put(-1,d)
        for i in xrange(n-1):
            L[i*n:(i+1)*n, i*n:(i+1)*n] = D
            L[(i+1)*n:(i+2)*n, i*n:(i+1)*n] = P
        # Last diagonal block
        L[-n:,-n:] = D
        return L

Here, ``put`` is sufficiently efficient that the benefit of constructing the
matrix by blocks is not apparent anymore. The slicing and block notation can
nevertheless be used for clarity. It could also be implemented as a combination
of ``find`` and ``put``, at the expense of memory consumption.

.. sourcecode:: ipython

   In [9]: %timeit -n10 -r3 L = poisson.poisson2d_sym_blk(1000)
   10 loops, best of 3: 246 ms per loop
   In [10]: %timeit -n10 -r3 L = poisson_vec.poisson2d_sym_blk_vec(1000)
   10 loops, best of 3: 232 ms per loop


Matlab Implementation
=====================

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

The Matlab functions above were place in a Matlab script names ``poisson.m``
which takes ``n`` as argument. It then calls ``poisson2d``, ``poisson2d_blk``
and ``poisson2d_kron`` successively, surrounding each call by ``tic`` and
``toc``. The tests were performed on a 2.4GHz Intel Core2 Duo running Matlab
7.6.0.324 (R2008a).

The results are as follows::

   >> poisson(100)
   poisson2d      Elapsed time is 1.731940 seconds.
   poisson2d_blk  Elapsed time is 0.804837 seconds.
   poisson2d_kron Elapsed time is 0.019118 seconds.

   >> poisson(300)
   poisson2d      Elapsed time is 145.979044 seconds.
   poisson2d_blk  Elapsed time is 32.785585 seconds.
   poisson2d_kron Elapsed time is 0.215165 seconds.

   >> poisson(500)
   poisson2d      Elapsed time is 2318.512099 seconds.
   poisson2d_blk  Elapsed time is 292.355093 seconds.
   poisson2d_kron Elapsed time is 0.627137 seconds.

   >> poisson(1000)
   poisson2d_kron Elapsed time is 2.317660 seconds.


It is striking to see how slow the straightforward ``poisson2d``
version is in Matlab. As we see in the next section, the Python version is
faster by several orders of magnitude.


Comparison with Matlab
======================

First, consider the simple ``Poisson2D`` function. The :ref:`table below <mpy1>`
summarizes the results of the previous section by giving timing ratios between
the Python and Matlab Poisson constructors.

.. _mpy1:
.. table:: Matlab vs. Python: Construction of 2D Poisson matrices.

   +-------+---------+--------+------------+---------------+-------------------+
   | ``n`` | Matlab  | Python | Python_vec | Matlab/Python | Matlab/Python_vec |
   +=======+=========+========+============+===============+===================+
   | 100   |    1.73 | 0.0382 | 0.00426    | 45.53         | 406.1             |
   +-------+---------+--------+------------+---------------+-------------------+
   | 300   |  145.98 | 0.3520 | 0.0317     | 414.72        | 4605.0            |
   +-------+---------+--------+------------+---------------+-------------------+
   | 500   | 2318.51 | 0.9800 | 0.0864     | 2365.8        | 26834.6           |
   +-------+---------+--------+------------+---------------+-------------------+
   | 1000  | |inf|   | 4.02   | 0.333      | |inf|         | |inf|             |
   +-------+---------+--------+------------+---------------+-------------------+

Unfortunately, since Matlab does not explicitly support symmetric matrices, we
cannot compare the other functions. For information only, we compare the block
version of the Python constructor with the Kronecker-product version of the
Matlab constructor. The results are in the :ref:`next table <mpy2>`.

.. _mpy2:
.. table:: Matlab vs. Python: Construction of 2D Poisson matrices---Fastest
           methods.

   +-------+---------+------------+-------------------+
   | ``n`` | Matlab  | Python_vec | Matlab/Python_vec |
   +=======+=========+============+===================+
   | 100   | 0.01912 | 0.0025     | 7.65              |
   +-------+---------+------------+-------------------+
   | 300   | 0.2152  | 0.0219     | 9.83              |
   +-------+---------+------------+-------------------+
   | 500   | 0.6271  | 0.0631     | 9.94              |
   +-------+---------+------------+-------------------+
   | 1000  | 2.318   | 0.232      | 9.99              |
   +-------+---------+------------+-------------------+

.. Shortcuts

.. |inf| replace:: :math:`\infty`
