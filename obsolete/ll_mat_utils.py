import re
import spmatrix

def ll_mat_from_mtx(f):
    """construct ll_mat matrix from MatrixMarket file

    parameters:

    f     file object, representing an opened file containing
          MatrixMarket data.
    """

    # parse header line
    line = f.readline()
    match = re.match('%%MatrixMarket\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)', line)
    if match == None:
        raise ValueError, "not a MatrixMarket file."
    if match.group(1) <> 'matrix' or match.group(2) <> 'coordinate' or match.group(3) <> 'real':
        raise NotImplementedError, 'matrix type not supported.'
    if match.group(4) == 'symmetric':
        isSymmetric = 1
    elif match.group(4) == 'unsymmetric':
        isSymmetric = 0
    else:
        raise ValueError, "invalid matrix type."
    # skip comment lines
    line = f.readline()
    while line[0] == '%':
        line = f.readline()
    # parse matrix dimension line
    match = re.match('\s*(\d+)\s+(\d+)\s+(\d+)', line)
    if match == None:
        raise ValueError, "invalid matrix size information in MatrixMarket file."
    dim = (int(match.group(1)), int(match.group(2)))
    nnz = int(match.group(3))
    if isSymmetric and dim[0] <> dim[1]:
        raise ValueError, 'symmteric matrix must be square'
    
    A = spmatrix.ll_mat(dim[0], dim[1])

    # read matrix elements
    k = 0
    while 1:
        lines = f.readlines(8*1024)
        if not lines: break
        for line in lines:
            elems = line.split()
            i = int(elems[0]) - 1
            j = int(elems[1]) - 1
            v = float(elems[2])
            A[i,j] = v
            if isSymmetric:
                A[j,i] = v
            k += 1
            if k == nnz:
                break
    f.close()
    return A
