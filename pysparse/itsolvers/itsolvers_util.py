"""This module defines encapsulation classes for all iterative solvers
implemented in the krylov module.

All classes provide a method "solve(b, x)" to approximatively solve the
linear system Ax=b using Krylov subspace methods.

The classes are intended to replace superlu.superlu_context in cases
where appropriate."""

from pysparse.sparse.pysparseMatrix import PysparseMatrix
from pysparse.itsolvers import krylov

class ItSolver:
    "Abstract base class for iteravtive solver classes."
    def __init__(self, A, debug):
        self.name = 'Generic'
        self.itsolver = None

        if isinstance(A, PysparseMatrix):
            self.A = A.matrix
        else:
            self.A = A

        self.nofCalled = 0
        self.totalIterations = 0
        self.lastIterations = 0
        self.lastInfo = 0
        self.debug = debug

    def set_name(self, name):
        self.name = name

    def set_solver(self, solver):
        self.itsolver = solver

    def solve(self, b, x, tol, maxit, K=None):
        "Solve Ax = b iteratively with zero initial guess"
        if not self.itsolver:
            raise NotImplementedError, 'This class cannot be used directly.'

        x[:] = 0
        info, iter, relres = self.itsolver(self.A, b, x, tol, maxit, K)
        self.nofCalled += 1
        self.totalIterations += iter
        self.lastIterations = iter
        self.lastInfo = info
        if self.debug:
            print 'iterative solver returned:', info, iter, relres
        if info < 0:
            raise RuntimeError, ('iterative solver %s returned error code %d' % (self.__class__.__name__, info), info, iter, relres)

    def __call__(self, *args, **kwargs):
        return self.solve(*args, **kwargs)

    def __str__(self):
        s = '<%s.%s instance>\n\n' % (self.__class__.__module__,
                                      self.__class__.__name__)
        attrs = ['nofCalled', 'totalIterations', 'lastIterations', 'lastInfo']
        for name in attrs:
            s += '   %s: %s\n' % (name, getattr(self, name))
        s += '\n'
        return s

    def __repr__(self):
        return self.__str__()


class Pcg(ItSolver):
    """
    Wrapper class for the preconditioned conjugate gradient iterative solver.
    """
    def __init__(self, A, debug=False):
        ItSolver.__init__(self, A, debug)
        self.set_name('PCG')
        self.set_solver(krylov.pcg)


class Minres(ItSolver):
    """
    Wrapper class for the minimum residual iterative solver.
    """
    def __init__(self, A, debug=False):
        ItSolver.__init__(self, A, debug)
        self.set_name('MINRES')
        self.set_solver(krylov.minres)


class Qmrs(ItSolver):
    """
    Wrapper class for the quasi-minimum residual smoothing iterative solver.
    """
    def __init__(self, A, debug=False):
        ItSolver.__init__(self, A, debug)
        self.set_name('QMRS')
        self.set_solver(krylov.qmrs)


class Cgs(ItSolver):
    """
    Wrapper class for the conjugate gradient squared iterative solver.
    """
    def __init__(self, A, debug=False):
        ItSolver.__init__(self, A, debug)
        self.set_name('CGS')
        self.set_solver(krylov.cgs)


class Bicgstab(ItSolver):
    """
    Wrapper class for the bi-conjugate gradient stabilized iterative solver.
    """
    def __init__(self, A, debug=False):
        ItSolver.__init__(self, A, debug)
        self.set_name('BICGSTAB')
        self.set_solver(krylov.bicgstab)


class Gmres(ItSolver):
    """
    Wrapper class for the generalized minimum residual iterative solver.
    """
    def __init__(self, A, debug=False):
        ItSolver.__init__(self, A, debug)
        self.set_name('GMRES')
        self.set_solver(krylov.gmres)


from pysparse.misc import Deprecated

@Deprecated('Use pysparse.itsolvers.Pcg instead.')
def pcg(*args, **kwargs):
    return  krylov.pcg(*args, **kwargs)

@Deprecated('Use pysparse.itsolvers.Minres instead.')
def minres(*args, **kwargs):
    return  krylov.minres(*args, **kwargs)

@Deprecated('Use pysparse.itsolvers.Qmrs instead.')
def qmrs(*args, **kwargs):
    return  krylov.qmrs(*args, **kwargs)

@Deprecated('Use pysparse.itsolvers.Cgs instead.')
def cgs(*args, **kwargs):
    return  krylov.cgs(*args, **kwargs)

@Deprecated('Use pysparse.itsolvers.Bicgstab instead.')
def gmres(*args, **kwargs):
    return  krylov.bicgstab(*args, **kwargs)

@Deprecated('Use pysparse.itsolvers.Gmres instead.')
def gmres(*args, **kwargs):
    return  krylov.gmres(*args, **kwargs)

if __name__ == '__main__':
    import numpy
    from pysparse.precon import precon
    from pysparse.tools import poisson

    n = 100
    n2 = n*n
    A = PysparseMatrix(matrix=poisson.poisson2d_sym(n))
    b = numpy.ones(n2); b /= numpy.linalg.norm(b)
    x = numpy.empty(n2)
    K = precon.ssor(A.matrix.to_sss())
    fmt = '%8s  %7.1e  %2d  %4d'

    def resid(A, b, x):
        r = b - A*x
        return numpy.linalg.norm(r)

    for Solver in [Pcg, Minres, Cgs, Qmrs, Gmres, Bicgstab]:
        solver = Solver(A)
        solver.solve(b, x, 1.0e-6, 3*n)
        print fmt % (solver.name, resid(A, b, x), solver.nofCalled,
                     solver.totalIterations)
        solver.solve(b, x, 1.0e-6, 3*n, K)
        print fmt % (solver.name, resid(A, b, x), solver.nofCalled,
                     solver.totalIterations)

