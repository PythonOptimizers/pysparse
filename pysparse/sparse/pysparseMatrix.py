#!/usr/bin/env python
"""
This module defines a few convenience classes as wrappers around ll_mat
objects. Being proper Python classes, they are subclassable. PysparseMatrix
objects have hooks for all methods of ll_mat objects.
"""

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "pysparseMatrix.py"
 #                                    created: 11/10/03 {3:15:38 PM}
 #                                last update: 1/3/07 {3:03:32 PM}
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #
 #  Description:
 #
 #  History
 #
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##

# A number of updates by Dominique Orban <dominique.orban@gmail.com>
# - allow creation of rectangular and square symmetric matrices
# - updates to __add__ and others to allow addition/subtraction of symmetric
#   matrices
# - new creator function PysparseMatrixSpDiags() to create banded matrices
#   with given diagonals.

__docformat__ = 'restructuredtext'

from pysparse.sparse import spmatrix
from pysparse.sparse.sparseMatrix import SparseMatrix
import numpy

class PysparseMatrix(SparseMatrix):
    """
    A PysparseMatrix is a class wrapper for the pysparse spmatrix sparse matrix
    type. This class facilitates matrix populating and allows intuitive
    operations on sparse matrices and vectors.

    :keywords:
        :nrow:       The number of rows of the matrix
        :ncol:       The number of columns of the matrix
        :size:       The common number of rows and columns, for a square matrix
        :bandwidth:  The bandwidth (if creating a band matrix)
        :matrix:     The starting `spmatrix` if there is one
        :sizeHint:   A guess on the number of nonzero elements of the matrix
        :symmetric:  A boolean indicating whether the matrix is symmetric.
        :storeZeros: A boolean indicating whether to store zero values.
    """

    def __init__(self, **kwargs):

        nrow = kwargs.get('nrow', 0)
        ncol = kwargs.get('ncol', 0)
        bandwidth = kwargs.get('bandwidth', 0)
        matrix = kwargs.get('matrix', None)
        sizeHint = kwargs.get('sizeHint', 0)
        storeZeros = kwargs.get('storeZeros', False)
        symmetric = 'symmetric' in kwargs and kwargs['symmetric']
        size = kwargs.get('size',0)

        if size > 0:
            if nrow > 0 or ncol > 0:
                if size != nrow or size != ncol:
                    msg =  'size argument was given but does not match '
                    msg += 'nrow and ncol'
                raise ValueError, msg
            else:
                nrow = ncol = size

        if matrix is not None:
            self.matrix = matrix
        else:
            if symmetric and nrow==ncol:
                if sizeHint is None:
                    sizeHint = nrow
                    if bandwidth > 0:
                        sizeHint += 2*(bandwidth-1)*(2*nrow-bandwidth-2)
                self.matrix = spmatrix.ll_mat_sym(nrow, sizeHint, storeZeros)
            else:
                if sizeHint is None:
                    sizeHint = min(nrow,ncol)
                    if bandwidth > 0:
                        sizeHint = bandwidth * (2*sizeHint-bandwidth-1)/2
                self.matrix = spmatrix.ll_mat(nrow, ncol, sizeHint, storeZeros)

    def isSymmetric(self):
        "Returns `True` is `self` is a symmetric matrix or `False` otherwise"
        if self.matrix.issym: return True
        return False

    def getNnz(self):
        "Returns the number of nonzero elements of `self`"
        return self.matrix.nnz

    def getMatrix(self):
        "Returns the underlying `ll_mat` sparse matrix of `self`"
        return self.matrix

    def copy(self):
        "Returns a (deep) copy of a sparse matrix"
        return PysparseMatrix(matrix = self.matrix.copy())

    def __coerce__(self, other):
        return self, other

    def __getattr__(self, name):
        if name == 'nnz':
            return self.getNnz()
        elif name == 'shape':
            return self.getShape()
        msg = 'No such attribute: %s' % name
        raise ValueError, msg

    def __getitem__(self, index):
        m = self.matrix[index]
        if isinstance(m, int) or isinstance(m, float):
            return m
        else:
            return PysparseMatrix(matrix = m, symmetric=self.matrix.issym)

    def __setitem__(self, index, value):
        #if type(value) is type(self):
        if isinstance(value, PysparseMatrix):
            self.matrix[index] = value.matrix
        else:
            self.matrix[index] = value

    def __iadd__(self, other):
        # In-place addition
        return self._iadd(self.getMatrix(), other)

    def _iadd(self, L, other, sign = 1):
        # In-place addition helper
        if self.isSymmetric() and not other.isSymmetric():
            L.generalize()
        L.shift(sign, other.getMatrix())
        return self

    def __add__(self, other):
        """
        Add two sparse matrices, return a new sparse matrix

            >>> L = PysparseMatrix(size = 3)
            >>> L.put([3.,10.,numpy.pi,2.5], [0,0,1,2], [2,1,1,0])
            >>> print L + PysparseIdentityMatrix(size = 3)
             1.000000  10.000000   3.000000
                ---     4.141593      ---
             2.500000      ---     1.000000

            >>> print L + 0
                ---    10.000000   3.000000
                ---     3.141593      ---
             2.500000      ---        ---

            >>> print L + 3
                ---    13.000000   6.000000
                ---     6.141593      ---
             5.500000      ---        ---
        """
        if other is 0 or other is 0.0:
            return self
        elif isinstance(other, int) or isinstance(other, float): #type(other) in [type(1), type(1.0)]:
            # Add give value to all elements of sparse matrix in nonzero pattern
            L = self.copy()
            val, irow, jcol = L.find()
            L.matrix.update_add_at( other*numpy.ones(val.shape), irow, jcol)
            return L
        elif type(self) is type(other):
            if self.getShape() != other.getShape():
                msg = 'Only sparse matrices of the same size may be added'
                raise TypeError, msg
            L = self.matrix.copy()
            if self.isSymmetric() and not other.isSymmetric():
                L.generalize()
            L.shift(1, other.getMatrix())
            return PysparseMatrix(matrix = L)

    def __sub__(self, other):

        if isinstance(other,int) or isinstance(other, float): #type(other) in [type(1), type(1.0)]:
            return self.__add__(-other)
        else:
            if self.getShape() != other.getShape():
                msg = 'Only sparse matrices of the same size may be subtracted'
                raise TypeError, msg
            L = self.matrix.copy()
            if self.isSymmetric() and not other.isSymmetric():
                L.generalize()
            L.shift(-1, other.getMatrix())
            return PysparseMatrix(matrix = L)

    def __isub__(self, other):
        # In-place subtraction
        return self._iadd(self.getMatrix(), other, sign=-1)

    def __mul__(self, other):
        """
        Multiply a sparse matrix by another sparse matrix

            >>> L1 = PysparseMatrix(size = 3)
            >>> L1.put([3.,10.,numpy.pi,2.5], [0,0,1,2], [2,1,1,0])
            >>> L2 = PysparseMatrix(size = 3)
            >>> L2.put(numpy.ones(3), numpy.arange(3), numpy.arange(3))
            >>> L2.put([4.38,12357.2,1.1], [2,1,0], [1,0,2])

            >>> tmp = numpy.array(((1.23572000e+05, 2.31400000e+01, 3.00000000e+00),
            ...                      (3.88212887e+04, 3.14159265e+00, 0.00000000e+00),
            ...                      (2.50000000e+00, 0.00000000e+00, 2.75000000e+00)))

            >>> numpy.allclose((L1 * L2).getNumpyArray(), tmp)
            1

        or a sparse matrix by a vector

            >>> tmp = numpy.array((29., 6.28318531, 2.5))
            >>> numpy.allclose(L1 * numpy.array((1,2,3),'d'), tmp)
            1

        or a vector by a sparse matrix

            >>> tmp = numpy.array((7.5, 16.28318531,  3.))
            >>> numpy.allclose(numpy.array((1,2,3),'d') * L1, tmp) ## The multiplication is broken. Numpy is calling __rmul__ for every element instead of with  the whole array.
            1
        """
        M, N = self.getShape()

        if isinstance(other, PysparseMatrix):
            if N != other.getShape()[0]:
                raise TypeError, 'Matrices dimensions do not match for product'

            p = spmatrix.matrixmultiply(self.matrix, other.getMatrix())
            return PysparseMatrix(matrix=p)
        else:
            shape = numpy.shape(other)
            if shape == ():  # other is a scalar
                p = self.matrix.copy()
                p.scale(other)
                return PysparseMatrix(matrix=p)
            elif shape == (N,):
                y = numpy.empty(M)
                self.matrix.matvec(other, y)
                return y
            else:
                raise TypeError, 'Cannot multiply objects'

    def __imul__(self, other):
        # In-place multiplication (by a scalar)
        #if type(other) not in [type(0), type(0.0)]:
        if not (isinstance(other, int) or isinstance(other, float)):
            raise TypeError, 'In-place multiplication is with scalars only'
        p = self.matrix
        p.scale(other)
        return self

    def __rmul__(self, other):
        # Compute  other * A  which is really  A^T * other
        if type(numpy.ones(1.0)) == type(other):
            M, N = self.getShape()
            y = numpy.empty(N)
            self.matrix.matvec_transp(other, y)
            return y
        else:
            return self * other

    def getShape(self):
        "Returns the shape ``(nrow,ncol)`` of a sparse matrix"
        return self.matrix.shape

    def col_scale(self, v):
        """
        Apply in-place column scaling. Each column is scaled by the
        corresponding component of ``v``, i.e., ``A[:,i] *= v[i]``.
        """
        return self.matrix.col_scale(v)

    def row_scale(self, v):
        """
        Apply in-place row scaling. Each row is scaled by the
        corresponding component of ``v``, i.e., ``A[i,:] *= v[i]``.

        """
        return self.matrix.row_scale(v)

    def find(self):
        """
        Returns three Numpy arrays to describe the sparsity pattern of ``self``
        in so-called coordinate (or triplet) format:

            >>> L = PysparseMatrix(size = 3)
            >>> L.put([3.,10.,numpy.pi,2.5], [0,0,1,2], [2,1,1,0])
            >>> (val,irow,jcol) = L.find()
            >>> val
            array([ 10.        ,   3.        ,   3.14159265,   2.5       ])
            >>> irow
            array([0, 0, 1, 2])
            >>> jcol
            array([1, 2, 1, 0])
        """
        return self.matrix.find()

    def put(self, value, id1=None, id2=None):
        """
        Put elements of ``value`` at positions of the matrix
        corresponding to ``(id1, id2)``

            >>> L = PysparseMatrix(size = 3)
            >>> L.put( [3.,10.,numpy.pi,2.5], [0,0,1,2], [2,1,1,0] )
            >>> print L
                ---    10.000000   3.000000
                ---     3.141593      ---
             2.500000      ---        ---
            >>> L.put(2*numpy.pi, range(3), range(3))
            >>> print L
             6.283185  10.000000   3.000000
                ---     6.283185      ---
             2.500000      ---     6.283185

        If ``value`` is a scalar, it has the same effect as the vector
        of appropriate length with all values equal to ``value``.
        If ``id1`` is omitted, it is replaced with ``range(nrow)``.
        If ``id2`` also is omitted, it is replaced with ``range(ncol)``.
        If ``id2`` is omitted but ``id1`` is present, ``id2`` is set to
        ``id1``.
        """
        nrow, ncol = self.getShape()
        if id2 is None and id1 is not None: id2 = id1
        if id1 is None:
            if not (isinstance(value, int) or isinstance(value, float)):
                id1 = numpy.arange(len(value), dtype=numpy.int)
            else:
                id1 = numpy.arange(nrow, dtype=numpy.int)
        if id2 is None:
            if not (isinstance(value, int) or isinstance(value, float)):
                id2 = numpy.arange(len(value), dtype=numpy.int)
            else:
                id2 = numpy.arange(ncol, dtype=numpy.int)
        self.matrix.put(value, id1, id2)
        return None

    def putDiagonal(self, vector):
        """
        Put elements of ``vector`` along diagonal of matrix

            >>> L = PysparseMatrix(size = 3)
            >>> L.putDiagonal([3.,10.,numpy.pi])
            >>> print L
             3.000000      ---        ---
                ---    10.000000      ---
                ---        ---     3.141593
            >>> L.putDiagonal([10.,3.])
            >>> print L
            10.000000      ---        ---
                ---     3.000000      ---
                ---        ---     3.141593
            >>> L.putDiagonal(2.7182)
            >>> print L
             2.718200      ---        ---
                ---     2.718200      ---
                ---        ---     2.718200

        """
        if isinstance(vector, int) or isinstance(vector, float):
            ids = numpy.arange(self.getShape()[0])
            #tmp = numpy.zeros((self.getShape()[0],), 'd')
            #tmp[:] = vector
            self.put(vector, ids, ids)
        else:
            #ids = numpy.arange(len(vector))
            self.matrix.put(vector) #, ids, ids)

    def take(self, id1=None, id2=None):
        """
        Extract elements at positions ``(irow[i], jcol[i])`` and place them in
        the array ``val``. In other words::

            for i in range(len(val)): val[i] = A[irow[i],jcol[i]]

        """
        nrow, ncol = self.getShape()
        if id2 is None and id1 is not None: id2 = id1
        if id1 is None:
            if not (isinstance(value, int) or isinstance(value, float)):
                id1 = numpy.arange(len(value), dtype=numpy.int)
            else:
                id1 = numpy.arange(nrow, dtype=numpy.int)
        if id2 is None:
            if not (isinstance(value, int) or isinstance(value, float)):
                id2 = numpy.arange(len(value), dtype=numpy.int)
            else:
                id2 = numpy.arange(ncol, dtype=numpy.int)
        vector = numpy.zeros(len(id1), 'd')
        self.matrix.take(vector, id1, id2)
        return vector

    def takeDiagonal(self):
        """
        Extract the diagonal of a matrix and place it in a Numpy array.
        """
        ids = numpy.arange(self.getShape()[0])
        return self.take(ids, ids)

    def addAt(self, vector, id1, id2):
        """
        Add elements of ``vector`` to the positions in the matrix corresponding
        to ``(id1,id2)``

            >>> L = PysparseMatrix(size = 3)
            >>> L.put([3.,10.,numpy.pi,2.5], [0,0,1,2], [2,1,1,0])
            >>> L.addAt((1.73,2.2,8.4,3.9,1.23), (1,2,0,0,1), (2,2,0,0,2))
            >>> print L
            12.300000  10.000000   3.000000
                ---     3.141593   2.960000
             2.500000      ---     2.200000
        """
        self.matrix.update_add_at(vector, id1, id2)

    def addAtDiagonal(self, vector):
        """
        Add the components of vector ``vector`` to the diagonal elements of the
        matrix.
        """
        #if type(vector) in [type(1), type(1.)]:
        if isinstance(vector, int) or isinstance(vector, float):
            ids = numpy.arange(self.getShape()[0])
            tmp = numpy.empty((self.getShape()[0],), 'd')
            tmp[:] = vector
            self.addAt(tmp, ids, ids)
        else:
            ids = numpy.arange(len(vector))
            self.addAt(vector, ids, ids)

    def getNumpyArray(self):
        """
        Convert a sparse matrix to a dense Numpy matrix.
        """
        shape = self.getShape()
        indices = numpy.indices(shape)
        numMatrix = self.take(indices[0].ravel(), indices[1].ravel())
        return numpy.reshape(numMatrix, shape)

    def matvec(self, x):
        """
        This method is required for scipy solvers.
        """
        return self * x

    def exportMmf(self, filename):
        """
        Exports the matrix to a Matrix Market file of the given filename.
        """
        self.matrix.export_mtx(filename)


class PysparseIdentityMatrix(PysparseMatrix):
    """
    Represents a sparse identity matrix for pysparse.

        >>> print PysparseIdentityMatrix(size = 3)
         1.000000      ---        ---
            ---     1.000000      ---
            ---        ---     1.000000
    """
    def __init__(self, size):

        PysparseMatrix.__init__(self, nrow=size, ncol=size,
                                bandwidth=1, symmetric=True)
        ids = numpy.arange(size)
        self.put(numpy.ones(size), ids, ids)


class PysparseSpDiagsMatrix(PysparseMatrix):
    """
    Represents a banded matrix with specified diagonals.

    *Example:* Create a tridiagonal matrix with 1's on the diagonal, 2's above
    the diagonal, and -2's below the diagonal.

        >>> from numpy import ones
        >>> e = ones(5)
        >>> print PysparseSpDiagsMatrix(size=5, vals=(-2*e,e,2*e), pos=(-1,0,1))
         1.000000   2.000000      ---        ---        ---
        -2.000000   1.000000   2.000000      ---        ---
            ---    -2.000000   1.000000   2.000000      ---
            ---        ---    -2.000000   1.000000   2.000000
            ---        ---        ---    -2.000000   1.000000

    Note that since the `pos[k]`-th diagonal has `size-|pos[k]|` elements, only
    that many first elements of `vals[k]` will be inserted.

    If the banded matrix is requested to be symmetric, elements above the
    main diagonal are not inserted.
    """
    def __init__(self, size, vals, pos, **kwargs):

        if type(pos) in [ type(()), type([]) ]: pos = numpy.array(pos)

        bw = max(numpy.abs(pos))
        diags = size - numpy.abs(pos)
        nz = sum(diags)
        kwargs.pop('bandwidth', True)
        kwargs.pop('sizeHint', True)
        PysparseMatrix.__init__(self, nrow=size, ncol=size,
                                bandwidth=bw, sizeHint=nz, **kwargs)

        # Insert elements on specified diagonals
        ndiags = len(pos)
        for k in range(ndiags):
            dk = diags[k]
            d = pos[k]
            if d >= 0 and not self.isSymmetric():
                self.put(vals[k][:dk], numpy.arange(dk), d + numpy.arange(dk))
            else:
                self.put(vals[k][:dk], -d + numpy.arange(dk), numpy.arange(dk))


# A subclass for use with Scipy solvers.
class PysparseMatrix4Scipy(PysparseMatrix):

    def __init__(self, **kwargs):
        PysparseMatrix.__init__(self, **kwargs)

    def matvec(self,x,y):
        """
        Return a new vector containing the matrix-vector product with `x`.
        This method is provided for compatibility with Scipy solvers.
        """
        return self.matrix.matvec(x,y)


def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
