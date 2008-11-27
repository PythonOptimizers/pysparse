.. Sparse matrix formats supplied by Pysparse
.. _formats-page:

=====================
Sparse Matrix Formats
=====================

This section describes the sparse matrix storage schemes available in
Pysparse. It also covers sparse matrix creation, population and conversion.

- Linked-list format (LL): a convenient format for creating and populating
  a sparse matrix, whether symmetric or general.
- Compressed sparse row format (CSR): a format designed to speed up
  matrix-vector products, but not well suited to matrix population and
  manipulation.
- Sparse Skyline format (SSS): a format for symmetric matrices designed to speed
  up matrix-vector products, but not well suited to matrix population and
  manipulation.

Linked-List Format
------------------

The linked-list format allows insertion and lookup of nonzero elements in
moderate time and without having to move too much data around. Internally, the
nonzero entries of a matrix are stored row by row in a linked list. Within
a given row, column indices are sorted in ascending order.

In Pysparse, matrices in linked-list format are created by using
the :class:`ll_mat` class.

This format resembles a sorted version of the coordinate format but with a data
structure that lends itself to fast insertion, removal and lookup.

Typically, a new matrix should be created as an :class:`ll_mat` and
populated. If necessary, it can then be converted to compressed sparse row or
sparse skyline format using the :meth:`to_csr` and :meth:`to_sss` methods.

The data structure for a matrix in linked-list format has the following
components:

``val``
    The double precision array ``val`` of length ``nalloc`` contains the
    non-zero entries of matrix.

``col``
    The integer array ``col`` of length ``nalloc`` contains the column indices
    of the non-zero entries stored in ``val``.

``link``
    the integer array ``link`` of length ``nalloc`` stores the pointer (index)
    to the next non-zero entry of the same row. A value of -1 indicates that
    there is no next entry.

``root``
    The integer array ``root`` of length ``n`` contains the pointers to the
    first entry of each row. The other entries of the same row can be located by
    following the ``link`` array.

``free``
    The integer ``free`` points to the first entry of the *free list*,
    i.e. a linked list of unoccupied spots in the ``val`` and ``col``
    arrays. This list is populated when non-zero entries are removed from the
    matrix.

Here ``n`` is the number of rows of the matrix and ``nalloc`` is number of
allocated elements in the arrays ``val``, ``col`` and ``link``. Note that the
number of nonzero entries stored is less than or equal to ``nalloc``, but the
``val``, ``col`` and ``link`` arrays can be enlarged dynamically if necessary.


Compressed Sparse Row Format
----------------------------

In CSR format, a sparse matrix is represented via three arrays:

``va``
    The double precision array ``va`` of length ``nnz`` contains the non-zero
    entries of the matrix, stored row by row.

``ja``
    The integer array ``ja`` of length ``nnz`` contains the column indices of
    the non-zero entries stored in ``va``.

``ia``
    The integer array ``ia`` of length ``n + 1`` contains the pointers (indices)
    to the beginning of each row in the arrays ``va`` and ``ja``. The last
    element of ``ia`` always has the value ``nnz + 1``.

Here ``n`` is the number of rows of the matrix and ``nnz`` is its number of
nonzero entries.

This format is particularly interesting for computing matrix-vector
products. Even though the order of the entries is not prescribed in this format,
we sort the entries of each row by ascending column indices. This enables us to
use more efficient algorithms for certain operations.


Sparse Skyline Format
---------------------

The SSS format is closely related to the CSR format. It is often used for sparse
*symmetric* matrices. The diagonal is stored in a separate (full) vector and the
strict lower triangle is stored in CSR format:

``va``
    The double precision array ``va`` of length ``nnz`` contains the non-zero
    entries of the strict lower triangle, stored row by row.

``ja``
    The integer array ``ja`` of length ``nnz`` contains the column indices of
    the non-zero entries stored in ``va``.

``ia``
    The integer array ``ia`` of length ``n + 1`` contains the pointers (indices)
    to the beginning of each row in the arrays ``va`` and ``ja``. The last
    element of ``ia`` always has the value ``nnz + 1``.

``da``
    The double precision array ``da`` of length ``n`` stores all diagonal
    entries of the matrix.

Here ``n`` is the order of the matrix and ``nnz`` is the number of nonzero
entries in the strict lower triangle.

We sort the entries of each row by ascending column indices, like we do with the
CSR format.  The SSS format has the advantage over the CSR format, that it
requires roughly half of the storage space and that the matrix-vector
multiplication can be implemented more efficiently
