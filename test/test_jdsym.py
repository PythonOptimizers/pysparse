import unittest
import numpy

from pysparse.itsolvers import krylov
from pysparse.eigen import jdsym
from pysparse.sparse import spmatrix

class JdsymTestCase(unittest.TestCase):

    def test_on_diagonal_matrix(self):
        a = spmatrix.ll_mat(3, 3)
        a.put([1, 2, 3])
        r = jdsym.jdsym(a, None, None, 3, 1, 1e-9, 100, krylov.qmrs)
        self.assertEqual(r[0], 3)
        self.assertTrue(numpy.allclose(r[1], [1, 2, 3]))
        self.assertTrue(numpy.allclose(numpy.abs(r[2]), numpy.identity(3)))

if __name__ == '__main__':
    unittest.main()
