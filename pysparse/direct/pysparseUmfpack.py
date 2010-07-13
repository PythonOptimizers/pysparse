"""
A framework for solving sparse linear systems of equations using an LU
factorization, by means of the unsymmetric multifrontal sparse LU factorization
package UMFPACK ([D04a]_, [D04b]_, [DD99]_, [DD97]_).

This package is appropriate for factorizing sparse square unsymmetric or
rectangular matrices.

See [UMF]_ for more information.

**References:**

.. [D04a] T. A. Davis, *A column pre-ordering strategy for the
          unsymmetric-pattern multifrontal method*, ACM Transactions on
          Mathematical Software, **30**\ (2), pp. 165-195, 2004.
.. [D04b] T. A. Davis, *Algorithm 832: UMFPACK, an unsymmetric-pattern
          multifrontal method*, ACM Transactions on Mathematical Software,
          **30**\ (2), pp. 196-199, 2004.
.. [DD99] T. A. Davis and I. S. Duff, *A combined unifrontal/multifrontal
          method for unsymmetric sparse matrices*, ACM Transactions on
          Mathematical Software, **25**\ (1), pp. 1-19, 1999.
.. [DD97] T. A. Davis and I. S. Duff, *An unsymmetric-pattern multifrontal
          method for sparse LU factorization*, SIAM Journal on Matrix Analysis
          and Applications, **18**\ (1), pp. 140-158, 1997.
.. [UMF] http://www.cise.ufl.edu/research/sparse/umfpack

"""

# To look into:
#  - wrap up other useful methods of UMFPACK
#  - allow other data types

__docformat__ = 'restructuredtext'

from pysparse.sparse import pysparseMatrix as psm
import numpy
#import resource

from pysparse.direct.directSolver import PysparseDirectSolver
from pysparse.direct import umfpack
from pysparse.tools.sptime import cputime
from string import upper

#def cputime():
#    return resource.getrusage(resource.RUSAGE_SELF)[0]

class PysparseUmfpackSolver( PysparseDirectSolver ):
    """
    `PysparseUmfpackSolver` is a wrapper class around the UMFPACK library for
    the factorization of full-rank n-by-m matrices. Only matrices with real
    coefficients are currently supported.

    :parameters:

        :A: A PysparseMatrix instance representing the matrix to be factorized.

    :keywords:

       :strategy: string that specifies what kind of ordering and pivoting
                  strategy UMFPACK should use. Valid values are 'auto',
                  'unsymmetric', 'symmetric' and '2by2'. Default: 'auto'

       :tol2by2: tolerance for the 2 by 2 strategy. Default: 0.1

       :scale: string that specifies the scaling UMFPACK should use. Valid
               values are 'none', 'sum', and 'max'. Default: 'sum'.

       :tolpivot: relative pivot tolerance for threshold partial pivoting with
                  row interchanges. Default: 0.1

       :tolsympivot: if diagonal pivoting is attempted, this parameter is used
                     to control when the diagonal is selected in a given pivot
                     column. Default: 0.0

    .. attribute:: LU

       An :class:`umfpack_context` object encapsulating the factorization.

    .. attribute:: sol

       The solution of the linear system after a call to :meth:`solve`.

    .. attribute:: L

       The L factor of the input matrix.

    .. attribute:: U

       The U factor of the input matrix.

    .. attribute:: P

       The row permutation used for the factorization.

    .. attribute:: Q

       The column permutation used for the factorization.

    .. attribute:: R

       The row scaling used during the factorization. See the documentation
       of :meth:`fetch_factors`.

    .. attribute:: factorizationTime

       The CPU time to perform the factorization.

    .. attribute:: solutionTime

       The CPU time to perform the forward and backward sweeps.

    .. attribute:: do_recip

       Nature of the row scaling. See :meth:`fetch_factors`.

    .. attribute:: lnz

       The number of nonzero elements in the factor L.

    .. attribute:: unz

       The number of nonzero elements in the factor U from which the diagonal
       was removed.

    .. attribute:: nz_udiag
 
       The number of nonzero elements on the diagonal of the factor U.
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
        self.lnz = self.unz = self.nz_udiag = None
        return

    def solve(self, rhs, **kwargs):
        """
        Solve the linear system  ``A x = rhs``. The result is
        placed in the :attr:`sol` member of the class instance.
        
        :parameters:
           :rhs: a Numpy vector of appropriate dimension.

        :keywords:
           :method: specifies the type of system being solved:
        
                    +-------------------+--------------------------------------+
                    |``"UMFPACK_A"``    | :math:`\mathbf{A} x = b` (default)   |
                    +-------------------+--------------------------------------+
                    |``"UMFPACK_At"``   | :math:`\mathbf{A}^T x = b`           |
                    +-------------------+--------------------------------------+
                    |``"UMFPACK_Pt_L"`` | :math:`\mathbf{P}^T \mathbf{L} x = b`|
                    +-------------------+--------------------------------------+
                    |``"UMFPACK_L"``    | :math:`\mathbf{L} x = b`             |
                    +-------------------+--------------------------------------+
                    |``"UMFPACK_Lt_P"`` | :math:`\mathbf{L}^T \mathbf{P} x = b`|
                    +-------------------+--------------------------------------+
                    |``"UMFPACK_Lt"``   | :math:`\mathbf{L}^T x = b`           |
                    +-------------------+--------------------------------------+
                    |``"UMFPACK_U_Qt"`` | :math:`\mathbf{U} \mathbf{Q}^T x = b`|
                    +-------------------+--------------------------------------+
                    |``"UMFPACK_U"``    | :math:`\mathbf{U} x = b`             |
                    +-------------------+--------------------------------------+
                    |``"UMFPACK_Q_Ut"`` | :math:`\mathbf{Q} \mathbf{U}^T x = b`|
                    +-------------------+--------------------------------------+
                    |``"UMFPACK_Ut"``   | :math:`\mathbf{U}^T x = b`           |
                    +-------------------+--------------------------------------+

           :irsteps: number of iterative refinement steps to attempt. Default: 2
        """
        method = kwargs.get('method', 'UMFPACK_A')
        irsteps = kwargs.get('irsteps', 2)
        if self.sol is None: self.sol = numpy.empty(self.ncol, self.type)
        t = cputime()
        self.LU.solve(rhs, self.sol, method, irsteps)
        self.solutionTime = cputime() - t
        return

    def fetch_lunz(self):
        """
        Retrieve the number of nonzeros in the factors. The results are stored
        in the members :attr:`lnz`, :attr:`unz` and :attr:`nz_udiag` of the
        class instance.
        """
        self.lnz, self.unz, self.nz_udiag = self.LU.lunz()

    def fetch_factors(self):
        """
        Retrieve the L and U factors of the input matrix along with the
        permutation matrices P and Q and the row scaling matrix R such that
 
        .. math:: \mathbf{P R A Q} = \mathbf{L U}.

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
