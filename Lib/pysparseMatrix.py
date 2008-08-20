#!/usr/bin/env python

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

from pysparse import spmatrix
from sparseMatrix import SparseMatrix
import numpy

class PysparseMatrix(SparseMatrix):
    
    """
    PysparseMatrix class wrapper for pysparse.
    PysparseMatrix is always NxN.
    Allows basic python operations __add__, __sub__ etc.
    Facilitate matrix populating in an easy way.

    """

    def __init__(self, **kwargs):
        """
        Creates a `PysparseMatrix`.

        :Currently accepted keywords include:

        - `nrow`: The number of rows of the matrix
        - `ncol`: The number of columns of the matrix
        - `bandwidth`: The bandwidth (if creating a band matrix)
        - `matrix`: The starting `spmatrix` if there is one
        - `sizeHint`: A guess on the number of nonzero elements of the matrix
        - `symmetric`: A boolean indicating whether the matrix is symmetric.
        """

        nrow = kwargs.get('nrow', 0)
        ncol = kwargs.get('ncol', 0)
        bandwidth = kwargs.get('bandwidth', 0)
        matrix = kwargs.get('matrix', None)
        sizeHint = kwargs.get('sizeHint', 0)
        symmetric = 'symmetric' in kwargs and kwargs['symmetric']

        if matrix is not None:
            self.matrix = matrix
        else:
            if symmetric and nrow==ncol:
                if sizeHint is None:
                    sizeHint = nrow
                    if bandwidth > 0:
                        sizeHint += 2*(bandwidth-1)*(2*nrow-bandwidth-2)
                self.matrix = spmatrix.ll_mat_sym(nrow, sizeHint)
            else:
                if sizeHint is None:
                    size = min(nrow,ncol)
                    sizeHint = size
                    if bandwidth > 0:
                        sizeHint = bandwidth * (2*size-bandwidth-1)/2
                self.matrix = spmatrix.ll_mat(nrow, ncol, sizeHint)

    def isSymmetric(self):
        if self.matrix.issym: return True
        return False

    def getNnz(self):
        return self.matrix.nnz

    def getMatrix(self):
        return self.matrix
    
    def copy(self):
        return PysparseMatrix(matrix = self.matrix.copy())
        
    def __coerce__(self, other): return self, other
    
    def __getattr__(self, name):
        if name == 'nnz':
            return self.getNnz()
        elif name == 'shape':
            return self.getShape()
        msg = 'No such attribute: %s' % name
        raise ValueError, msg

    def __getitem__(self, index):
        m = self.matrix[index]
        if type(m) is type(0) or type(m) is type(0.):
            return m
        else:
            return PysparseMatrix(matrix = m, symmetric=self.matrix.issym)

    def __iadd__(self, other):
        # To implement  L += K
        return self._iadd(self.getMatrix(), other)
        
    def _iadd(self, L, other, sign = 1):
        if other != 0:
            if self.isSymmetric() and not other.isSymmetric():
                L.generalize()
            L.shift(sign, other.getMatrix())
        return self

    def __add__(self, other):
        """
        Add two sparse matrices
        
            >>> L = PysparseMatrix(size = 3)
            >>> L.put((3.,10.,numpy.pi,2.5), (0,0,1,2), (2,1,1,0))
            >>> print L + PysparseIdentityMatrix(size = 3)
             1.000000  10.000000   3.000000  
                ---     4.141593      ---    
             2.500000      ---     1.000000  
             
            >>> print L + 0
                ---    10.000000   3.000000  
                ---     3.141593      ---    
             2.500000      ---        ---    
            
            >>> print L + 3
            Traceback (most recent call last):
            ...
            AttributeError: 'int' object has no attribute 'getMatrix'
        """

        if self.getShape() != other.getShape():
            raise TypeError, 'Only sparse matrices of same size may be added'
        if other is 0:
            return self
        else:
            L = self.matrix.copy()
            if self.isSymmetric() and not other.isSymmetric():
                L.generalize()
            L.shift(1, other.getMatrix())
            return PysparseMatrix(matrix = L)
        
    def __sub__(self, other):

        if self.getShape() != other.getShape():
            raise TypeError, 'Only sparse matrices of same size may be subtracted'

        if other is 0:
            return self
        else:
            L = self.matrix.copy()
            if self.isSymmetric() and not other.isSymmetric():
                L.generalize()
            L.shift(-1, other.getMatrix())
            return PysparseMatrix(matrix = L)

    def __isub__(self, other):
        # To implement L -= K
        return self._iadd(self.getMatrix(), other, -1)

    def __mul__(self, other):
        """
        Multiply a sparse matrix by another sparse matrix
        
            >>> L1 = PysparseMatrix(size = 3)
            >>> L1.put((3.,10.,numpy.pi,2.5), (0,0,1,2), (2,1,1,0))
            >>> L2 = PysparseIdentityMatrix(size = 3)
            >>> L2.put((4.38,12357.2,1.1), (2,1,0), (1,0,2))
            
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
        N = self.matrix.shape[1]

        if isinstance(other, PysparseMatrix):
            if N != other.getShape()[0]:
                raise TypeError, 'Matrices dimensions do not match for product'

            return PysparseMatrix(matrix = spmatrix.matrixmultiply(self.matrix, other.getMatrix()))
        else:
            shape = numpy.shape(other)
            if shape == ():
                L = spmatrix.ll_mat(N, N, N)
                L.put(other * numpy.ones(N))
                return PysparseMatrix(matrix = spmatrix.matrixmultiply(self.matrix, L))
            elif shape == (N,):
                y = numpy.empty(N) #other.copy()
                self.matrix.matvec(other, y)
                return y
            else:
                raise TypeError, 'Cannot multiply objects'
            
    def __rmul__(self, other):
        if type(numpy.ones(1)) == type(other):
            y = numpy.empty(numpy.shape(other)) #other.copy()
            self.matrix.matvec_transp(other, y)
            return y
        else:
            return self * other
            
    def getShape(self):
        return self.matrix.shape
        
    def put(self, vector, id1, id2):
        """
        Put elements of `vector` at positions of the matrix corresponding to (`id1`, `id2`)
        
            >>> L = PysparseMatrix(size = 3)
            >>> L.put((3.,10.,numpy.pi,2.5), (0,0,1,2), (2,1,1,0))
            >>> print L
                ---    10.000000   3.000000  
                ---     3.141593      ---    
             2.500000      ---        ---    
        """
        self.matrix.put(vector, id1, id2)

    def putDiagonal(self, vector):
        """
        Put elements of `vector` along diagonal of matrix
        
            >>> L = PysparseMatrix(size = 3)
            >>> L.putDiagonal((3.,10.,numpy.pi))
            >>> print L
             3.000000      ---        ---    
                ---    10.000000      ---    
                ---        ---     3.141593  
            >>> L.putDiagonal((10.,3.))
            >>> print L
            10.000000      ---        ---    
                ---     3.000000      ---    
                ---        ---     3.141593  
        """
        if type(vector) in [type(1), type(1.)]:
            ids = numpy.arange(self.getShape()[0])
            tmp = numpy.zeros((self.getShape()[0],), 'd')
            tmp[:] = vector
            self.put(tmp, ids, ids)
        else:
            ids = numpy.arange(len(vector))
            self.put(vector, ids, ids)

    def take(self, id1, id2):
        vector = numpy.zeros(len(id1), 'd')
        self.matrix.take(vector, id1, id2)
        return vector

    def takeDiagonal(self):
        ids = numpy.arange(self.getShape()[0])
        return self.take(ids, ids)

    def addAt(self, vector, id1, id2):
        """
        Add elements of `vector` to the positions in the matrix corresponding to (`id1`,`id2`)
        
            >>> L = PysparseMatrix(size = 3)
            >>> L.put((3.,10.,numpy.pi,2.5), (0,0,1,2), (2,1,1,0))
            >>> L.addAt((1.73,2.2,8.4,3.9,1.23), (1,2,0,0,1), (2,2,0,0,2))
            >>> print L
            12.300000  10.000000   3.000000  
                ---     3.141593   2.960000  
             2.500000      ---     2.200000  
        """
        self.matrix.update_add_at(vector, id1, id2)

    def addAtDiagonal(self, vector):
        if type(vector) in [type(1), type(1.)]:
            ids = numpy.arange(self.getShape()[0])
            tmp = numpy.empty((self.getShape()[0],), 'd')
            tmp[:] = vector
            self.addAt(tmp, ids, ids)
        else:
            ids = numpy.arange(len(vector))
            self.addAt(vector, ids, ids)

    def getNumpyArray(self):
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
    """
    def __init__(self, size):
        """
        Create a sparse matrix with '1' in the diagonal
        
            >>> print PysparseIdentityMatrix(size = 3)
             1.000000      ---        ---    
                ---     1.000000      ---    
                ---        ---     1.000000  
        """
        PysparseMatrix.__init__(self, nrow=size, ncol=size,
                                bandwidth=1, symmetric=True)
        ids = numpy.arange(size)
        self.put(numpy.ones(size), ids, ids)

class PysparseSpDiagsMatrix(PysparseMatrix):
    """
    Represents a banded matrix with specified diagonals.
    """
    def __init__(self, size, vals, pos, **kwargs):
        """
        Example:
        Create a tridiagonal matrix with 1's on the diagonal, 2's above the
        diagonal, and -2's below the diagonal.

            >>> from numpy import ones
            >>> e = ones(5)
            >>> print PysparseSpDiagsMatrix(size=5, vals=(-2*e,e,2*e), pos=(-1,0,1))

        Note that since the pos[k]-th diagonal has size-|pos[k]| elements, only
        that many first elements of vals[k] will be inserted.

        If the banded matrix is requested to be symmetric, elements above the
        main diagonal are not inserted.
        """
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

        
def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
