"""
A framework for solving sparse linear systems of equations using an LU
factorization, by means of the supernodal sparse LU factorization package
SuperLU.

This package is appropriate for factorizing sparse square unsymmetric or
rectangular matrices.

See http://crd.lbl.gov/~xiaoye/SuperLU for more information.

References:

- J. W. Demmel, S. C. Eisenstat, J. R. Gilbert, X. S. Li and J. W. H. Liu,
  *A supernodal approach to sparse partial pivoting*, SIAM Journal on Matrix
  Analysis and Applications **20**(3), pp. 720-755, 1999.
- J. W. Demmel, J. R. Gilbert and X. S. Li, *An Asynchronous Parallel Supernodal
  Algorithm for Sparse Gaussian Elimination*, SIAM Journal on Matrix Analysis
  and Applications **20**(4), pp. 915-952, 1999.
- X. S. Li and J. W. Demmel, *SuperLU_DIST: A Scalable Distributed-Memory Sparse
  Direct Solver for Unsymmetric Linear Systems*, ACM Transactions on
  Mathematical Software **29**(2), pp. 110-140, 2003.
"""

# To look into:
#  - allow other data types
#  - rely on user-installed SuperLU library?

__docformat__ = 'restructuredtext'

import pysparseMatrix as psm
import numpy
import resource

from directSolver import PysparseDirectSolver
from pysparse  import superlu

def cputime():
    return resource.getrusage(resource.RUSAGE_SELF)[0]

class PysparseSuperLUSolver( PysparseDirectSolver ):
    """
    `PysparseSuperLUSolver` is a wrapper class around the SuperLu library for
    the factorization of full-rank n-by-m matrices. Only matrices with real
    coefficients are currently supported.

    The input matrix A should be supplied as a PysparseMatrix instance.

    Currently accepted keywords include

    `symmetric`
      a boolean indicating that the user wishes to use symmetric mode. In
      symmetric mode, permc_spec=2 must be chosen and diag_pivot_thresh must be
      small, e.g., 0.0 or 0.1. Since the value of diag_pivot_thresh is up to the
      user, setting symmetric to True does *not* automatically set permc_spec
      and diag_pivot_thresh to appropriate values.

    `diag_pivot_thresh`
      a float value between 0 and 1 representing the threshold for partial
      pivoting (0 = no pivoting, 1 = always perform partial pivoting).
      Default: 1.0.

    `drop_tol`
      the value of a drop tolerance, between 0 and 1, if an incomplete
      factorization is desired (0 = exact factorization). This keyword has no
      effect if using SuperLU version 2.0 and below.
      Default: 0.0

    `relax`
      an integer controling the degree of relaxing supernodes.
      Default: 1.

    `panel_size`
      an integer specifying the maximum number of columns to form a panel.
      Default: 10.

    `permc_spec`
      an integer specifying the ordering strategy used during the factorization.
      0 : natural ordering,
      1 : MMD applied to the structure of A^T * A
      2 : MMD applied to the structure of A^T + A
      3 : COLAMD.
      Default: 2.

    """
    def __init__(self, A, **kwargs):
        PysparseDirectSolver.__init__(self, A, **kwargs)

        self.type = numpy.float
        self.nrow, self.ncol = A.getShape()
        t = cputime()
        self.LU = superlu.factorize(A.matrix.to_csr(), **kwargs)
        self.factorizationTime = cputime() - t
        self.solutionTime = 0.0
        self.sol = None
        self.L = self.U = None
        return

    def solve(self, rhs, transpose = False):
        """
        `solve(rhs)`:
        Solve the linear system  ``A x = rhs``, where ``A`` is the input matrix
        and ``rhs`` is a Numpy vector of appropriate dimension. The result is
        placed in the ``sol`` member of the class instance.

        If the optional argument ``transpose`` is ``True``, the transpose system
        ``A^T x = rhs`` is solved.
        """
        if self.sol is None: self.sol = numpy.empty(self.ncol, self.type)
        transp = 'N'
        if transpose: transp = 'T'
        t = cputime()
        self.LU.solve(rhs, self.sol, transp)
        self.solutionTime = cputime() - t
        return

    def fetch_lunz(self):
        """
        `fetch_lunz()`:
        Retrieve the number of nonzeros in the factors L and U together. The
        result is stored in the member ``lunz`` of the class instance.
        """
        self.lunz = self.LU.nnz

    def fetch_factors(self):
        """
        Not yet available.
        """
        raise NotImplementedError
