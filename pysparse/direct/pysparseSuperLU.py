"""
A framework for solving sparse linear systems of equations using an LU
factorization, by means of the supernodal sparse LU factorization package
SuperLU ([DEGLL99]_, [DGL99]_, [LD03]_).

This package is appropriate for factorizing sparse square unsymmetric or
rectangular matrices.

See [SLU]_ for more information.

**References:**

.. [DEGLL99] J. W. Demmel, S. C. Eisenstat, J. R. Gilbert, X. S. Li and
             J. W. H. Liu, *A supernodal approach to sparse partial pivoting*,
             SIAM Journal on Matrix Analysis and Applications **20**\ (3),
             pp. 720-755, 1999.
.. [DGL99] J. W. Demmel, J. R. Gilbert and X. S. Li,
           *An Asynchronous Parallel Supernodal Algorithm for Sparse Gaussian
           Elimination*, SIAM Journal on Matrix Analysis and Applications
           **20**\ (4), pp. 915-952, 1999.
.. [LD03] X. S. Li and J. W. Demmel, *SuperLU_DIST: A Scalable
          Distributed-Memory Sparse Direct Solver for Unsymmetric Linear
          Systems*, ACM Transactions on Mathematical Software **29**\ (2),
          pp. 110-140, 2003.
.. [SLU] http://crd.lbl.gov/~xiaoye/SuperLU

"""

# To look into:
#  - allow other data types

__docformat__ = 'restructuredtext'

from pysparse.sparse import pysparseMatrix as psm
import numpy
#import resource

from pysparse.direct.directSolver import PysparseDirectSolver
from pysparse.direct import superlu
from pysparse.tools.sptime import cputime

#def cputime():
#    return resource.getrusage(resource.RUSAGE_SELF)[0]

class PysparseSuperLUSolver( PysparseDirectSolver ):
    """
    `PysparseSuperLUSolver` is a wrapper class around the SuperLu library for
    the factorization of full-rank n-by-m matrices. Only matrices with real
    coefficients are currently supported.

    :parameters:

       :A: The matrix to be factorized, supplied as a PysparseMatrix instance.

    :keywords:

       :symmetric: a boolean indicating that the user wishes to use symmetric
                   mode. In symmetric mode, ``permc_spec=2`` must be chosen and
                   ``diag_pivot_thresh`` must be small, e.g., 0.0 or 0.1. Since
                   the value of ``diag_pivot_thresh`` is up to the user, setting
                   ``symmetric`` to ``True`` does *not* automatically set
                   ``permc_spec`` and ``diag_pivot_thresh`` to appropriate
                   values.

       :diag_pivot_thresh: a float value between 0 and 1 representing the
                           threshold for partial pivoting (0 = no pivoting,
                           1 = always perform partial pivoting). Default: 1.0.

       :drop_tol: the value of a drop tolerance, between 0 and 1, if an
                  incomplete factorization is desired (0 = exact factorization).
                  This keyword does not exist if using SuperLU version 2.0 and
                  below. In more recent version of SuperLU, the keyword is
                  accepted but has no effect. Default: 0.0

       :relax: an integer controling the degree of relaxing supernodes.
               Default: 1.

       :panel_size: an integer specifying the maximum number of columns to form
                    a panel. Default: 10.

       :permc_spec: an integer specifying the ordering strategy used during the
                    factorization.

                    0. natural ordering,
                    1. MMD applied to the structure of
                       :math:`\mathbf{A}^T \mathbf{A}`
                    2. MMD applied to the structure of
                       :math:`\mathbf{A}^T + \mathbf{A}`
                    3. COLAMD.

                    Default: 2.

    .. attribute:: LU

       A :class:`superlu_context` object encapsulating the factorization.

    .. attribute:: sol

       The solution of the linear system after a call to :meth:`solve`.

    .. attribute:: factorizationTime

       The CPU time to perform the factorization.

    .. attribute:: solutionTime

       The CPU time to perform the forward and backward sweeps.

    .. attribute:: lunz

       The number of nonzero elements in the factors L and U together after a
       call to :meth:`fetch_lunz`.
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
        Solve the linear system  ``A x = rhs``, where ``A`` is the input matrix
        and ``rhs`` is a Numpy vector of appropriate dimension. The result is
        placed in the :attr:`sol` member of the class instance.

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
        Retrieve the number of nonzeros in the factors L and U together. The
        result is stored in the member :attr:`lunz` of the class instance.
        """
        self.lunz = self.LU.nnz

    def fetch_factors(self):
        """
        Not yet available.
        """
        raise NotImplementedError
