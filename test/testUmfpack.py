import math, unittest
import numpy
from pysparse.tools import poisson
from pysparse.sparse.pysparseMatrix import PysparseMatrix, PysparseIdentityMatrix, PysparseSpDiagsMatrix
from pysparse.direct.pysparseUmfpack import PysparseUmfpackSolver

import sys

def macheps():
    "compute machine epsilon"
    eps = 1.0
    while (1.0 + eps > 1.0):
        eps /= 2.0
    return 2.0 * eps
    
class SpDiagsTestCase(unittest.TestCase):
    def setUp(self):
        self.n = 10000
        self.e = numpy.ones(self.n)
        self.A = PysparseSpDiagsMatrix(size=self.n,
                                       vals=(-2*self.e,self.e,2*self.e),
                                       pos=(-1,0,1))
        self.b = self.A * self.e
        self.tol = 100*macheps()
        self.LU = None
        self.err = 0.0
        self.fmt = '\t%8.2e  %8.2e  %8d  %8d  %8d  %6.2f  %6.2f\n'

    def tearDown(self):
        self.LU.fetch_lunz()
        sys.stdout.write(self.fmt % (self.err, self.tol, self.A.getNnz(),
                                     self.LU.lnz, self.LU.unz,
                                     self.LU.factorizationTime,
                                     self.LU.solutionTime))
        
    def computeError(self, x):
        self.err = numpy.linalg.norm(x - self.e, ord=numpy.inf)
        return self.err

    def testVanilla(self):
        self.LU = PysparseUmfpackSolver(self.A)
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testUnsymmetric(self):
        self.LU = PysparseUmfpackSolver(self.A, strategy='unsymmetric')
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testSymmetric(self):
        self.LU = PysparseUmfpackSolver(self.A, strategy='symmetric')
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def test2by2(self):
        self.LU = PysparseUmfpackSolver(self.A, strategy='2by2')
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testNoScaling(self):
        self.LU = PysparseUmfpackSolver(self.A, scale='none')
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testScaleMax(self):
        self.LU = PysparseUmfpackSolver(self.A, scale='max')
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)


            
class Poisson1dTestCase(unittest.TestCase):
    def setUp(self):
        self.n = 50000
        self.B = PysparseMatrix( matrix=poisson.poisson1d(self.n) )
        
        self.x_exact = numpy.ones(self.n)/math.sqrt(self.n)
        self.normx = 1.0/math.sqrt(self.n)
        self.b = self.B * self.x_exact
        
        lmbd_min = 4.0 * math.sin(math.pi/2.0/self.n) ** 2
        lmbd_max = 4.0 * math.sin((self.n - 1)*math.pi/2.0/self.n) ** 2
        cond = lmbd_max/lmbd_min
        self.tol = cond * macheps()
        self.relerr = 0.0
        self.nnz = self.B.getNnz()
        self.LU = None
        self.fmt = '\t%8.2e  %8.2e  %8d  %8d  %8d  %6.2f  %6.2f\n'

    def computeError(self, x):
        absErr = numpy.linalg.norm(x - self.x_exact, ord=numpy.inf)
        self.relerr = absErr/(1 + self.normx)
        return self.relerr

    def tearDown(self):
        self.LU.fetch_lunz()
        sys.stdout.write(self.fmt % (self.relerr, self.tol, self.nnz,
                                     self.LU.lnz, self.LU.unz,
                                     self.LU.factorizationTime,
                                     self.LU.solutionTime))
        
    def testPoisson1dDefault(self):
        self.LU = PysparseUmfpackSolver(self.B)
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testPoisson1dUnsymmetric(self):
        self.LU = PysparseUmfpackSolver(self.B, strategy='unsymmetric')
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testPoisson1dSymmetric(self):
        self.LU = PysparseUmfpackSolver(self.B, strategy='symmetric')
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testPoisson1d2by2(self):
        self.LU = PysparseUmfpackSolver(self.B, strategy='2by2')
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testPoisson1dNoScaling(self):
        self.LU = PysparseUmfpackSolver(self.B, scale='none')
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testPoisson1dScaleMax(self):
        self.LU = PysparseUmfpackSolver(self.B, scale='max')
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)



        
class Poisson2dTestCase(unittest.TestCase):
    def setUp(self):
        self.n = 200
        self.B = PysparseMatrix( matrix=poisson.poisson2d(self.n) )
        
        self.x_exact = numpy.ones(self.n*self.n)/self.n
        self.normx = 1.0/self.n
        self.b = self.B * self.x_exact

        h = 1.0/self.n
        lmbd_min = 4.0/h/h * (math.sin(math.pi*h/2.0) ** 2 +
                              math.sin(math.pi*h/2.0) ** 2)
        lmbd_max = 4.0/h/h * (math.sin((self.n - 1)*math.pi*h/2.0) ** 2 +
                              math.sin((self.n - 1)*math.pi*h/2.0) ** 2)
        cond = lmbd_max/lmbd_min
        self.tol = cond * macheps()
        self.relerr = 0.0
        self.nnz = self.B.getNnz()
        self.LU = None
        self.fmt = '\t%8.2e  %8.2e  %8d  %8d  %8d  %6.2f  %6.2f\n'

    def computeError(self, x):
        absErr = numpy.linalg.norm(x - self.x_exact, ord=numpy.inf)
        self.relerr = absErr/(1 + self.normx)
        return self.relerr

    def tearDown(self):
        self.LU.fetch_lunz()
        sys.stdout.write(self.fmt % (self.relerr, self.tol, self.nnz,
                                     self.LU.lnz, self.LU.unz,
                                     self.LU.factorizationTime,
                                     self.LU.solutionTime))
        
    def testPoisson2dDefault(self):
        self.LU = PysparseUmfpackSolver(self.B)
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testPoisson2dUnsymmetric(self):
        self.LU = PysparseUmfpackSolver(self.B, strategy='unsymmetric')
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testPoisson2dSymmetric(self):
        self.LU = PysparseUmfpackSolver(self.B, strategy='symmetric')
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testPoisson2d2by2(self):
        self.LU = PysparseUmfpackSolver(self.B, strategy='2by2')
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testPoisson2dNoScaling(self):
        self.LU = PysparseUmfpackSolver(self.B, scale='none')
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testPoisson2dScaleMax(self):
        self.LU = PysparseUmfpackSolver(self.B, scale='max')
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)



if __name__ == '__main__':
    headfmt = '\t%8s  %8s  %8s  %8s  %8s  %6s  %6s'
    header = headfmt % ('RelErr','Tol','nnz(A)','nnz(L)',
                        'nnz(U)','Fact','Solve')
    lhead = len(header)
    sys.stderr.write(header + '\n')
    sys.stderr.write('\t' + '-' * lhead + '\n')
    unittest.main()
