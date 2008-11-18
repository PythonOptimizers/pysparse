.. Description of the higher-level PysparseMatrix class
.. versionadded:: 1.0.1
.. _pysparsematrix-page:

==================================
Higher-Level Sparse Matrix Classes
==================================

The :mod:`pysparseMatrix` module
--------------------------------

.. automodule:: pysparseMatrix

.. autoclass:: PysparseMatrix
   :show-inheritance: 
   :members: 
   :inherited-members: 
   :undoc-members:

Creating an Identity Matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: PysparseIdentityMatrix
   :show-inheritance: 
   :members: 
   :undoc-members:

Creating Sparse Matrices from Diagonals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: PysparseSpDiagsMatrix
   :show-inheritance: 
   :members: 
   :undoc-members: 

Fancy Indexing
^^^^^^^^^^^^^^

Fancy indexing carries over to ``PysparseMatrix`` objects and is used exactly in
the same way as with ``ll_mat`` objects. Refer to Section :ref:`spmatrix-page`
for details.
