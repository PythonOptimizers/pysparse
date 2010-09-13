from pysparse.sparse import spmatrix
from pysparse.direct import umfpack
import numpy as np
from poisson2dvec import poisson2d_vec_sym_blk

l = spmatrix.ll_mat(5, 5)
l[0,0] = 2.0
l[0,1] = 3.0
l[1,4] = 6.0
l[1,0] = 3.0
l[1,2] = 4.0
l[2,1] = -1.0
l[2,2] = -3.0
l[2,3] = 2.0
l[3,2] = 1.0
l[4,1] = 4.0
l[4,2] = 2.0
l[4,4] = 1.0

b = np.array([8.0, 45.0, -3.0, 3.0, 19.0], "d")
x = np.zeros(5, "d")
umf = umfpack.factorize(l)
umf.solve(b, x, 'UMFPACK_A')
print umf.getlists()
print x

print "------------------------------"

n = 50
L = poisson2d_vec_sym_blk(n)
b = np.ones(n * n, 'd')
x = np.empty(n * n, 'd')
umf = umfpack.factorize(L)
umf.solve(b, x, 'UMFPACK_A')

r = np.empty(n * n, 'd')
L.matvec(x, r)
r = b - r
print 'norm(b - A * x) = %f' % np.linalg.norm(r)
