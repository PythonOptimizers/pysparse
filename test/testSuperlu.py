import math, unittest
import numpy
from pysparse.tools import poisson
from pysparse.sparse.pysparseMatrix import PysparseMatrix, PysparseIdentityMatrix, PysparseSpDiagsMatrix
from pysparse.direct.pysparseSuperLU import PysparseSuperLUSolver

import sys

def macheps():
    "compute machine epsilon"
    eps = 1.0
    while (1.0 + eps > 1.0):
        eps /= 2.0
    return 2.0 * eps

#
# Tests that are commented out work fine but take a long time...
#
    
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
        self.relerr = 0.0
        self.descr = ''
        self.fmt = '\t%10s  %8.2e  %8.2e  %8d  %8d  %6.2f  %6.2f\n'

    def tearDown(self):
        self.LU.fetch_lunz()
        sys.stdout.write(self.fmt % (self.descr, self.relerr, self.tol,
                                     self.A.getNnz(),
                                     self.LU.lunz, self.LU.factorizationTime,
                                     self.LU.solutionTime))
        
    def computeError(self, x):
        self.relerr = numpy.linalg.norm(x - self.e, ord=numpy.inf)
        return self.relerr

    #def testVanilla(self):
    #    self.descr = 'spdgs-dflt'
    #    self.LU = PysparseSuperLUSolver(self.A)
    #    self.LU.solve(self.b)
    #    self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testTrivialThresh(self):
        self.descr = 'spdgs-trsh'
        self.LU = PysparseSuperLUSolver(self.A, diag_pivot_thresh=0.5)
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    #def testTrivialRelax(self):
    #    self.descr = 'spdgs-relx'
    #    self.LU = PysparseSuperLUSolver(self.A, relax=20)
    #    self.LU.solve(self.b)
    #    self.failUnless(self.computeError(self.LU.sol) < self.tol)

    #def testTrivialPanel(self):
    #    self.descr = 'spdgs-size'
    #    self.LU = PysparseSuperLUSolver(self.A, panel_size=1)
    #    self.LU.solve(self.b)
    #    self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testTrivialperm0(self):
        self.descr = 'spdgs-prm0'
        self.LU = PysparseSuperLUSolver(self.A, permc_spec=0)
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testTrivialperm1(self):
        self.descr = 'spdgs-prm1'
        self.LU = PysparseSuperLUSolver(self.A, permc_spec=1)
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    #def testTrivialperm2(self):
    #    self.descr = 'spdgs-prm2'
    #    self.LU = PysparseSuperLUSolver(self.A, permc_spec=2)
    #    self.LU.solve(self.b)
    #    self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testTrivialperm3(self):
        self.descr = 'spdgs-prm3'
        self.LU = PysparseSuperLUSolver(self.A, permc_spec=3)
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)
            
class Poisson1dTestCase(unittest.TestCase):

    def setUp(self):
        self.n = 50000
        self.A = PysparseMatrix(matrix=poisson.poisson1d_sym(self.n))
        
        self.x_exact = numpy.ones(self.n)/math.sqrt(self.n)
        self.normx = 1.0/math.sqrt(self.n)
        self.b = self.A * self.x_exact
        
        lmbd_min = 4.0 * math.sin(math.pi/2.0/self.n) ** 2
        lmbd_max = 4.0 * math.sin((self.n - 1)*math.pi/2.0/self.n) ** 2
        cond = lmbd_max/lmbd_min
        self.tol = cond * macheps()
        self.LU = None
        self.relerr = 0.0
        self.descr = ''
        self.fmt = '\t%10s  %8.2e  %8.2e  %8d  %8d  %6.2f  %6.2f\n'

    def tearDown(self):
        self.LU.fetch_lunz()
        sys.stdout.write(self.fmt % (self.descr, self.relerr, self.tol,
                                     self.A.getNnz(),
                                     self.LU.lunz, self.LU.factorizationTime,
                                     self.LU.solutionTime))
        
    def computeError(self, x):
        absErr = numpy.linalg.norm(x-self.x_exact, ord=numpy.inf)
        self.relerr = absErr/(1+self.normx)
        return self.relerr

    def testPoisson1dDefault(self):
        self.descr = 'poi1d-dflt'
        self.LU = PysparseSuperLUSolver(self.A)
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testPoisson1dThresh(self):
        self.descr = 'poi1d-trsh'
        self.LU = PysparseSuperLUSolver(self.A, diag_pivot_thresh=0.5)
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testPoisson1dRelax(self):
        self.descr = 'poi1d-relx'
        self.LU = PysparseSuperLUSolver(self.A, relax=20)
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testPoisson1dPanel(self):
        self.descr = 'poi1d-size'
        self.LU = PysparseSuperLUSolver(self.A, panel_size=1)
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testPoisson1dperm0(self):
        self.descr = 'poi1d-prm0'
        self.LU = PysparseSuperLUSolver(self.A, permc_spec=0)
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testPoisson1dperm1(self):
        self.descr = 'poi1d-prm1'
        self.LU = PysparseSuperLUSolver(self.A, permc_spec=1)
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testPoisson1dperm2(self):
        self.descr = 'poi1d-prm2'
        self.LU = PysparseSuperLUSolver(self.A, permc_spec=2)
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testPoisson1dperm3(self):
        self.descr = 'poi1d-prm3'
        self.LU = PysparseSuperLUSolver(self.A, permc_spec=3)
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)
        
class Poisson2dTestCase(unittest.TestCase):
    def setUp(self):
        self.n = 200
        self.A = PysparseMatrix(matrix=poisson.poisson2d_sym_blk(self.n))
        
        self.x_exact = numpy.ones(self.n*self.n)/self.n
        self.normx = 1.0/self.n
        self.b = self.A * self.x_exact

        h = 1.0 / self.n
        lmbd_min = 4.0/h/h * (math.sin(math.pi*h/2.0) ** 2 +
                              math.sin(math.pi*h/2.0) ** 2)
        lmbd_max = 4.0/h/h * (math.sin((self.n - 1)*math.pi*h/2.0) ** 2 +
                              math.sin((self.n - 1)*math.pi*h/2.0) ** 2)
        cond = lmbd_max/lmbd_min
        self.tol = cond * macheps()
        self.LU = None
        self.relerr = 0.0
        self.descr = ''
        self.fmt = '\t%10s  %8.2e  %8.2e  %8d  %8d  %6.2f  %6.2f\n'

    def tearDown(self):
        self.LU.fetch_lunz()
        sys.stdout.write(self.fmt % (self.descr, self.relerr, self.tol,
                                     self.A.getNnz(),
                                     self.LU.lunz, self.LU.factorizationTime,
                                     self.LU.solutionTime))
        
    def computeError(self, x):
        absErr = numpy.linalg.norm(x-self.x_exact, ord=numpy.inf)
        self.relerr = absErr/(1+self.normx)
        return self.relerr

    def testPoisson2dDefault(self):
        self.descr = 'poi2d-dftl'
        self.LU = PysparseSuperLUSolver(self.A)
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testPoisson2dThresh(self):
        self.descr = 'poi2d-trsh'
        self.LU = PysparseSuperLUSolver(self.A, diag_pivot_thresh=0.5)
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testPoisson2dRelax(self):
        self.descr = 'poi2d-relx'
        self.LU = PysparseSuperLUSolver(self.A, relax=20)
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testPoisson2dPanel(self):
        self.descr = 'poi2d-size'
        self.LU = PysparseSuperLUSolver(self.A, panel_size=1)
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testPoisson2dperm0(self):
        self.descr = 'poi2d-prm0'
        self.LU = PysparseSuperLUSolver(self.A, permc_spec=0)
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testPoisson2dperm1(self):
        self.descr = 'poi2d-prm1'
        self.LU = PysparseSuperLUSolver(self.A, permc_spec=1)
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testPoisson2dperm2(self):
        self.descr = 'poi2d-prm2'
        self.LU = PysparseSuperLUSolver(self.A, permc_spec=2)
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)

    def testPoisson2dperm3(self):
        self.descr = 'poi2d-prm3'
        self.LU = PysparseSuperLUSolver(self.A, permc_spec=3)
        self.LU.solve(self.b)
        self.failUnless(self.computeError(self.LU.sol) < self.tol)
        
if __name__ == '__main__':
    headfmt = '\t%10s  %8s  %8s  %8s  %8s  %6s  %6s'
    header = headfmt % ('Test','RelErr','Tol','nnz(A)',
                        'nnz(L+U)','Fact','Solve')
    lhead = len(header)
    sys.stderr.write(header + '\n')
    sys.stderr.write('\t' + '-' * lhead + '\n')
    unittest.main()
