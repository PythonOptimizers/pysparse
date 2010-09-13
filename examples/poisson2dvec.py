# Poisson 2D constructors. Illustrate vectorization.
from pysparse.sparse import spmatrix
import numpy as np

def poisson2d_vec(n):
    n2 = n*n
    L = spmatrix.ll_mat(n2, n2, 5*n2-4*n)
    d = np.arange(n2, dtype=np.int)
    L.put(4.0, d)
    L.put(-1.0, d[:-n], d[n:])
    L.put(-1.0, d[n:], d[:-n])
    for i in xrange(n):
        di = d[i*n:(i+1)*n]
        L.put(-1.0, di[1:], di[:-1])
        L.put(-1.0, di[:-1], di[1:])
    return L

def poisson2d_vec_sym(n):
    n2 = n*n
    L = spmatrix.ll_mat_sym(n2, 3*n2-2*n)
    d = np.arange(n2, dtype=np.int)
    L.put(4.0, d)
    L.put(-1.0, d[n:], d[:-n])
    for i in xrange(n):
        di = d[i*n:(i+1)*n]
        L.put(-1.0, di[:-1], di[1:])
    return L

def poisson2d_vec_sym_blk(n):
    n2 = n*n
    L = spmatrix.ll_mat_sym(n2, 3*n2-2*n)
    D = spmatrix.ll_mat_sym(n, 2*n-1)
    d = np.arange(n, dtype=np.int)
    D.put(4.0, d)
    D.put(-1.0, d[1:], d[:-1])
    P = spmatrix.ll_mat_sym(n, n-1)
    P.put(-1,d)
    for i in xrange(n-1):
        L[i*n:(i+1)*n, i*n:(i+1)*n] = D
        L[(i+1)*n:(i+2)*n, i*n:(i+1)*n] = P
    # Last diagonal block
    L[-n:,-n:] = D
    return L
