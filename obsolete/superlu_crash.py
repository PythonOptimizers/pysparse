import superlu, spmatrix, Numeric

A = spmatrix.ll_mat_from_mtx("/home/geus/matrices/superlu_crash_2x2.mtx")
A = A.to_csr()
A
LU = superlu.factorize(A)

x = Numeric.zeros(2, 'd')
b = Numeric.ones(2, 'd')

LU.solve(b, x)
print x

A.matvec(x,b)
print b
