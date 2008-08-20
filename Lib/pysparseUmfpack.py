"""
A framework for solving sparse linear systems of equations using an LU
factorization, by means of the unsymmetric multifrontal sparse LU factorization
package UMFPACK.

This package is appropriate for factorizing sparse square unsymmetric or
rectangular matrices.

See http://www.cise.ufl.edu/research/sparse/umfpack for more information.

References:

- T. A. Davis, *A column pre-ordering strategy for the unsymmetric-pattern
  multifrontal method*, ACM Transactions on Mathematical Software, **30**(2),
  pp. 165-195, 2004.
- T. A. Davis, *Algorithm 832: UMFPACK, an unsymmetric-pattern multifrontal
  method*, ACM Transactions on Mathematical Software, **30**(2), pp. 196-199,
  2004.
- T. A. Davis and I. S. Duff, *A combined unifrontal/multifrontal method for
  unsymmetric sparse matrices*, ACM Transactions on Mathematical Software,
  **25**(1), pp. 1-19, 1999.
- T. A. Davis and I. S. Duff, *An unsymmetric-pattern multifrontal method for
  sparse LU factorization*, SIAM Journal on Matrix Analysis and Applications,
  **18**(1), pp. 140-158, 1997.
"""

# To look into:
#  - wrap up other useful methods of UMFPACK
#  - allow other data types
#  - rely on user-installed UMFPACK library?

__docformat__ = 'restructuredtext'

import pysparseMatrix as psm
import numpy
import resource

from directSolver import PysparseDirectSolver
from pysparse     import umfpack
from string       import upper

def cputime():
    return resource.getrusage(resource.RUSAGE_SELF)[0]

class PysparseUmfpackSolver( PysparseDirectSolver ):
    """
    `PysparseUmfpackSolver` is a wrapper class around the UMFPACK library for
    the factorization of full-rank n-by-m matrices. Only matrices with real
    coefficients are currently supported.

    The input matrix A should be supplied as a PysparseMatrix instance.

    Currently accepted keywords include

    `strategy`
      string that specifies what kind of ordering and pivoting strategy UMFPACK
      should use. Valid values are 'auto', 'unsymmetric', 'symmetric' and
      '2by2'. Default: 'auto'

    `tol2by2`
      tolerance for the 2 by 2 strategy. Default: 0.1

    `scale`
      string that specifies the scaling UMFPACK should use. Valid values are
      'none', 'sum', and 'max'. Default: 'sum'.

    `tolpivot`
      relative pivot tolerance for threshold partial pivoting with row
      interchanges. Default: 0.1

    `tolsympivot`
      if diagonal pivoting is attempted, this parameter is used to control when
      the diagonal is selected in a given pivot column. Default: 0.0

    `irstep`
      number of iterative refinement steps to attempt. Default: 2
    """
    def __init__(self, A, **kwargs):
        PysparseDirectSolver.__init__(self, A, **kwargs)

        if 'strategy' in kwargs.keys():
            strategy = upper(kwargs.get('strategy'))
            if strategy not in ['AUTO', 'UNSYMMETRIC', 'SYMMETRIC', '2BY2']:
                strategy = 'AUTO'
            kwargs['strategy'] = 'UMFPACK_STRATEGY_' + strategy

        if 'scale' in kwargs.keys():
            scale = upper(kwargs.get('scale'))
            if scale not in ['NONE', 'SUM', 'MAX']: scale = 'SUM'
            kwargs['scale'] = 'UMFPACK_SCALE_' + scale
        
        self.type = numpy.float
        self.nrow, self.ncol = A.getShape()
        t = cputime()
        self.LU = umfpack.factorize(A.matrix, **kwargs)
        self.factorizationTime = cputime() - t
        self.solutionTime = 0.0
        self.sol = None
        self.L = self.U = None
        self.P = self.Q = self.R = None
        self.do_recip = False
        return

    def solve(self, rhs, method='UMFPACK_A'):
        """
        `solve(rhs)`:
        Solve the linear system  ``A x = rhs``, where ``A`` is the input matrix
        and ``rhs`` is a Numpy vector of appropriate dimension. The result is
        placed in the ``sol`` member of the class instance.

        The optional ``method`` argument specifies the type of system being
        solved:
        
        - ``"UMFPACK_A"``    : Solve A x = b     (default)
        - ``"UMFPACK_At"``   : Solve A^t x = b
        - ``"UMFPACK_Pt_L"`` : Solve P^T L x = b
        - ``"UMFPACK_L"``    : Solve L x = b
        - ``"UMFPACK_Lt_P"`` : Solve L^t P x = b
        - ``"UMFPACK_Lt"``   : Solve L^t x = b
        - ``"UMFPACK_U_Qt"`` : Solve U Q^t x = b
        - ``"UMFPACK_U"``    : Solve U x = b
        - ``"UMFPACK_Q_Ut"`` : Solve Q U^t x = b
        - ``"UMFPACK_Ut"``   : Solve U^t x = b

        """
        if self.sol is None: self.sol = numpy.empty(self.ncol, self.type)
        t = cputime()
        self.LU.solve(rhs, self.sol, method)
        self.solutionTime = cputime() - t
        return

    def fetch_lunz(self):
        """
        `fetch_lunz()`:
        Retrieve the number of nonzeros in the factors. The results are stored
        in the members ``lnz``, ``unz`` and ``nz_udiag`` of the class instance.
        """
        self.lnz, self.unz, self.nz_udiag = self.LU.lunz()

    def fetch_factors(self):
        """
        `fetch_factors()`:
        Retrieve the L and U factors of the input matrix along with the
        permutation matrices P and Q and the row scaling matrix R such that
          P R A Q = L U.

        The matrices P, R and Q are stored as Numpy arrays. L and U are stored
        as PysparseMatrix instances and are lower triangular and upper
        triangular, respectively.

        R is a row-scaling diagonal matrix such that

        - the i-th row of A has been divided    by R[i] if ``do_recip = True``,
        - the i-th row of A has been multiplied by R[i] if ``do_recip = False``.

        """
        (L, U, self.P, self.Q, self.R, self.do_recip) = self.LU.lu()
        self.L = psm.PysparseMatrix(matrix=L)
        self.U = psm.PysparseMatrix(matrix=U)
        return
