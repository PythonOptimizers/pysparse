# Benchmark the pcg module of PySparse implementing
# a preconditioned conjugate gradient.
# Compare different preconditioners.

from pysparse.sparse import spmatrix
from pysparse.itsolvers.krylov import pcg
from pysparse.precon import precon
import numpy as np
import resource
import sys
import os

def usage():
    progname = sys.argv[0]
    sys.stderr.write('Usage\n');
    sys.stderr.write('  %-s problem [problem ...]\n' % progname)
    sys.stderr.write('  where each problem is the file name of a matrix\n')
    sys.stderr.write('  in MatrixMarket sparse format (*.mtx).\n')

def cputime():
        return resource.getrusage(resource.RUSAGE_SELF)[0]

def test_pcg(ProblemList, tol=1.0e-6):

    if len(ProblemList) == 0:
        usage()
        sys.exit(1)
    
    header1 = '%10s  %6s  %6s  ' % ('Name', 'n', 'nnz')
    header2 = '%6s  %8s  %8s  %4s  %6s  %6s\n' % ('iter','relres','error','info','form M','solve')
    dheader1 = '%10s  %6d  %6d  '
    dheader2 = '%6d  %8.1e  %8.1e  %4d  %6.2f  %6.2f\n' 
    lhead1 = len(header1)
    lhead2 = len(header2)
    lhead = lhead1 + lhead2
    sys.stderr.write('-' * lhead + '\n')
    sys.stderr.write(header1)
    sys.stderr.write(header2)
    sys.stderr.write('-' * lhead + '\n')

    # Record timings for each preconditioner
    timings = { 'None' : [],
                'Diagonal' : [],
                'SSOR' : []
              }

    for problem in ProblemList:
    
        A = spmatrix.ll_mat_from_mtx(problem)
        (m, n) = A.shape
        if m != n: break
    
        prob = os.path.basename(problem)
        if prob[-4:] == '.mtx': prob = prob[:-4]

        # Right-hand side is Ae
        e = np.ones(n, 'd')
        b = np.empty(n, 'd')
        A.matvec(e, b)

        sys.stdout.write(dheader1 % (prob, n, A.nnz))

        # No preconditioner
        x = np.zeros(n, 'd')
        t = cputime()
        info, iter, relres = pcg(A, b, x, tol, 2*n)
        t_noprec = cputime() - t
        err = np.linalg.norm(x-e, ord=np.Inf)
        sys.stdout.write(dheader2 % (iter, relres, err, info, 0.0, t_noprec))
        timings['None'].append(t_noprec)

        # Diagonal preconditioner
        x = np.zeros(n, 'd')
        t = cputime()
        M = precon.jacobi(A, 1.0, 1)
        t_getM_diag = cputime() - t
        t = cputime()
        info, iter, relres = pcg(A, b, x, tol, 2*n, M)
        t_diag = cputime() - t
        err = np.linalg.norm(x-e, ord=np.Inf)
        sys.stdout.write(lhead1 * ' ')
        sys.stdout.write(dheader2 % (iter,relres,err,info,t_getM_diag,t_diag))
        timings['Diagonal'].append(t_diag)

        # SSOR preconditioner
        # It appears that steps=1 and omega=1.0 are nearly optimal in all cases
        x = np.zeros(n, 'd')
        t = cputime()
        M = precon.ssor(A.to_sss(), 1.0, 1)
        t_getM_ssor = cputime() - t
        t = cputime()
        info, iter, relres = pcg(A, b, x, tol, 2*n, M)
        t_ssor = cputime() - t
        err = np.linalg.norm(x-e, ord=np.Inf)
        sys.stdout.write(lhead1 * ' ')
        sys.stdout.write(dheader2 % (iter,relres,err,info,t_getM_ssor,t_ssor))
        timings['SSOR'].append(t_ssor)
        sys.stderr.write('-' * lhead + '\n')

    return timings


if __name__ == '__main__':
    problem_list = sys.argv[1:]
    timings = test_pcg(problem_list)

    try:
        import matplotlib
        if matplotlib.__version__ < '0.65':
            import matplotlib.matlab as MM
        else:
            import matplotlib.pylab as MM
    except:
        print ' If you had Matplotlib installed, you would be looking'
        print ' at timing plots right now...'
        sys.exit(0)
        
    darkblue = '#2c11cf'
    lightblue = '#8f84e0'
    steelblue = '#5d82ef'
    colors = [ darkblue, lightblue, steelblue ]

    style0 = '-k'
    style1 = '--k'
    style2 = '-.k'
    styles = [ style0, style1, style2 ]

    # For the number of iterations, use first value of p as reference
    x = range(len(problem_list))
    ax = MM.subplot(111)
    lgnd = []
    i = 0
    for k in timings.keys():
        ax.plot(x, timings[k], color = colors[i], linewidth = 2)
        lgnd.append(k)
        i += 1
        if i > len(colors): i = 0
    ax.legend(lgnd, 'upper left')
    ax.set_title('Solve time for each preconditioner type in PCG')
    #ax.set_xticklabels(problem_list, rotation = 45, horizontalalignment = 'right')
    MM.show()
