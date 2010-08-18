"""
A framework for solving sparse linear systems of equations using a direct
factorization.
"""
__docformat__ = 'restructuredtext'

from pysparse.sparse import pysparseMatrix as psm
import numpy

class PysparseDirectSolver:
    """
    `PysparseDirectSolver` is a generic class and should be subclassed.
    """
    def __init__(self, A, **kwargs):
        return

    def solve(self, b, **kwargs):
        raise NotImplementedError
        return

