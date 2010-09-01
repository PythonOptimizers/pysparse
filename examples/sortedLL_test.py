from pysparse.sparse.spmatrix import ll_mat
import numpy as np
import time

n = 1000
nnz = 50000
A = ll_mat(n, n, nnz)

R = np.random.random_integers(0, n-1, (nnz,2))

t1 = time.clock()

for k in xrange(nnz):
    A[int(R[k,0]), int(R[k,1])] = k
    
print 'Time for populating matrix: %8.2f s' % (time.clock() - t1)

print 'nnz(A) = ', A.nnz

B = A[:,:]
A.shift(-1.0, B)
print A

