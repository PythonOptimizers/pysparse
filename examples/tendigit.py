"""
Solves problem 7 of the One Hundred Dollars, One Hundred Digits Challenge.
"""

import numpy as np
from pysparse.sparse import spmatrix
from pysparse.itsolvers.krylov import minres
from pysparse.precon import precon

def get_primes(nofPrimes):
    primes = np.zeros(nofPrimes, 'i')
    primes[0] = 2
    nof = 1
    i = 3
    while 1:
        for p in primes[:nof]:
            if i%p == 0 or p*p > i: break
        if i%p <> 0:
            primes[nof] = i
            nof += 1
            if nof >= nofPrimes:
                break
        i = i+2
    return primes

n = 20000
print 'Generating first %d primes...' % n
primes = get_primes(n)

print 'Assembling coefficient matrix...'
A = spmatrix.ll_mat_sym(n, n*8)
d = 1
while d < n:
    for i in range(d, n):
        A[i,i-d] = 1.0
    d *= 2
for i in range(n):
    A[i,i] = 1.0 * primes[i]

A = A.to_sss()
K = precon.ssor(A)

print 'Solving linear system...'
b = np.zeros(n); b[0] = 1.0
x = np.empty(n)
info, iter, relres = minres(A, b, x, 1e-16, n, K)

print info, iter, relres
print '%.16e' % x[0]
