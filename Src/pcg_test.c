/* PCG_TEST - test code for PCG code
 */

#include <assert.h>
#include <blas.h>
#include <fortran.h>
#include <lapack.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pcg.h"

/* matrix A */
static int n_s;
static double *va_s, *da_s;
static int *ja_s, *ia_s;

/* diag preconditioner */
static double *idiag_s;

static void randvec(int n, int k, double *x)
{
  int IDIST, m;
  static int ISEED[4] = {2, 3, 5, 7};
  int i;

  for (i = 0; i < k; i ++) {
    /* call fortran lapack routine dlanrv to a get random vector */
    IDIST = 1;
    m = n;
    F77(dlarnv)(&IDIST, ISEED, &m, x+i*n);    
  }
}

void read_msr(char *fname, int *n,
              double **va, double **da, int **ja, int **ia) {
  int i, j, k, nnz;
  FILE *f;

  f = fopen(fname, "r");
  assert(f != NULL);
  fscanf(f, "%d", n);

  (*da) = (double *)malloc(*n * sizeof(double));
  (*ia) = (int *)malloc((*n+1) * sizeof(int));
  assert(da && ia);

  for (i = 0; i <= *n; i ++) {
    fscanf(f, "%d", &((*ia)[i]));
    (*ia)[i] --;
  }
  nnz = (*ia)[*n];

  (*va) = (double *)malloc((nnz * sizeof(double)));
  (*ja) = (int *)malloc(nnz * sizeof(int));
  assert(va && ja);

  for (k = 0; k < nnz; k ++) {
    fscanf(f, "%d", &j);
    (*ja)[k] = j-1;
  }
  for (k = 0; k < nnz; k ++)
    fscanf(f, "%lf", &((*va)[k]));

  for (i = 0; i < *n; i ++)
    fscanf(f, "%lf", &((*da)[i]));

  fclose(f);
}

/* MATVEC - matrix vector multiplications
 */
void matveca(double *x, double *y) {
  int i, k;

  for(i = 0; i < n_s; i ++) {
    y[i] = da_s[i]*x[i];

    for(k = ia_s[i]; k < ia_s[i+1]; k ++) {
      y[i] = y[i] + va_s[k]*x[ja_s[k]];
      y[ja_s[k]] = y[ja_s[k]] + va_s[k]*x[i];
    }
  }
}

/* PREC_DIAG - diagonal preconditioner
 */
void prec_diag(double *x, double *y) {
  int i;
  for (i = 0; i < n_s; i ++)
    y[i] = idiag_s[i]*x[i];
}


void main(void) {
  double relres, nrmx;
  int n, i, iter, flag, ONE=1;
  double *x, *x1, *b, *work;

  read_msr("/home/geus/jdbsym/test/edge6x3x5_B.msr", 
	   &n_s, &va_s, &da_s, &ja_s, &ia_s);
  n = n_s;

  x       = (double *)malloc(n * sizeof(double));
  x1      = (double *)malloc(n * sizeof(double));
  b       = (double *)malloc(n * sizeof(double));
  work    = (double *)malloc(4*n * sizeof(double));
  idiag_s = (double *)malloc(n * sizeof(double));
  assert(x && x1 && b && work && idiag_s);
  
  randvec(n, 1, b);
  for (i = 0; i < n; i ++)
    b[i] = 1.0;
  matveca(b, x1);
  
  for (i = 0; i < n; i ++)
    x[i] = 0.0;

  for (i = 0; i < n; i ++)
    idiag_s[i] = 1.0 / da_s[i];
  
  pcg(n, 
      x, 
      b,
      1e-11, 
      500,
      1,
      &iter, 
      &relres, 
      &flag,
      work,
      matveca,
      prec_diag);
  
  printf("iter = %d, relres = %e, flag = %d\n",
	 iter, relres, flag);

  printf("Should converge in 161 iterations.\n");

  nrmx = F77(dnrm2)(&n, b, &ONE);
  printf("norm(b) = %e\n", nrmx);
  nrmx = F77(dnrm2)(&n, x, &ONE);
  printf("norm(x) = %e\n", nrmx);
  matveca(x, work);
  for (i = 0; i < n; i ++)
    work[i] = b[i] - work[i];
  nrmx = F77(dnrm2)(&n, work, &ONE);
  printf("norm(b - A*x) = %e\n", nrmx);
  nrmx = F77(dnrm2)(&n, work, &ONE) / F77(dnrm2)(&n, b, &ONE);
  printf("norm(b - A*x)/norm(b) = %e\n", nrmx);
}
