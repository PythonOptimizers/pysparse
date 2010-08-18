# A vectorized Poisson module
from pysparse.sparse import spmatrix
import numpy

def poisson1d_vec(n):
    L = spmatrix.ll_mat(n, n, 3*n-2)
    e = numpy.ones(n)
    d = numpy.arange(n, dtype=numpy.int)
    L.put(2*e, d, d)
    L.put(-e[1:], d[1:], d[:-1])
    L.put(-e[1:], d[:-1], d[1:])
    return L

def poisson1d_sym_vec(n):
    L = spmatrix.ll_mat_sym(n, 2*n-1)
    e = numpy.ones(n)
    d = numpy.arange(n, dtype=numpy.int)
    L.put(2*e, d, d)
    L.put(-e[1:], d[1:], d[:-1])
    return L

def poisson2d_vec(n):
    # First version, proceed block by block
    n2 = n*n
    L = spmatrix.ll_mat(n2, n2, 5*n2-4*n)
    e = numpy.ones(n)
    d = numpy.arange(n, dtype=numpy.int)
    din = d
    for i in xrange(n):
        # Diagonal blocks
        L.put(4*e, din, din)
        L.put(-e[1:], din[1:], din[:-1])
        L.put(-e[1:], din[:-1], din[1:])
        # Outer blocks
        L.put(-e, n+din, din)
        L.put(-e, din, n+din)
        din = d + i*n
    # Last diagonal block
    L.put(4*e, din, din)
    L.put(-e[1:], din[1:], din[:-1])
    L.put(-e[1:], din[:-1], din[1:])
    return L

def poisson2d_vec2(n):
    # Second version, allocate long arrays
    n2 = n*n
    L = spmatrix.ll_mat(n2, n2, 5*n2-4*n)
    e = numpy.ones(n2)
    d = numpy.arange(n2, dtype=numpy.int)
    L.put(4*e, d, d)
    din = d[:n]
    for i in xrange(n):
        # Diagonal blocks
        L.put(-e[:n-1], din[1:], din[:-1])
        L.put(-e[:n-1], din[:-1], din[1:])
        # Outer blocks
        L.put(-e[:n], n+din, din)
        L.put(-e[:n], din, n+din)
        din = d[i*n:(i+1)*n]
    # Last diagonal block
    L.put(-e[:n-1], din[1:], din[:-1])
    L.put(-e[:n-1], din[:-1], din[1:])
    return L

def poisson2d_sym_vec(n):
    n2 = n*n
    L = spmatrix.ll_mat_sym(n2, 3*n2-2*n)
    e = numpy.ones(n)
    d = numpy.arange(n, dtype=numpy.int)
    din = d
    for i in xrange(n):
        # Diagonal blocks
        L.put(4*e, din, din)
        L.put(-e[1:], din[1:], din[:-1])
        # Outer blocks
        L.put(-e, n+din, din)
        din = d + i*n
    # Last diagonal block
    L.put(4*e, din, din)
    L.put(-e[1:], din[1:], din[:-1])
    return L

def poisson2d_sym_blk_vec(n):
    n2 = n*n
    L = spmatrix.ll_mat_sym(n2, 3*n2-2*n)
    D = spmatrix.ll_mat_sym(n, 2*n-1)
    e = numpy.ones(n)
    d = numpy.arange(n, dtype=numpy.int)
    D.put(4*e, d, d)
    D.put(-e[1:], d[1:], d[:-1])
    P = spmatrix.ll_mat(n, n, n-1)
    P.put(-e,d,d)
    for i in xrange(n-1):
        L[i*n:(i+1)*n, i*n:(i+1)*n] = D
        L[(i+1)*n:(i+2)*n, i*n:(i+1)*n] = P
    # Last diagonal block
    L[n2-n:n2, n2-n:n2] = D
    return L
