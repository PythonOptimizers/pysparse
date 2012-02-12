/*
This module provides an interface to the UMFPACK library. The
license of the UMFPACK library is available below.

UMFPACK Version 4.1 (Apr. 30, 2003),  Copyright (c) 2003 by Timothy A.
Davis.  All Rights Reserved.

UMFPACK License:

    Your use or distribution of UMFPACK or any modified version of
    UMFPACK implies that you agree to this License.

    THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
    EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

    Permission is hereby granted to use or copy this program, provided
    that the Copyright, this License, and the Availability of the original
    version is retained on all copies.  User documentation of any code that
    uses UMFPACK or any modified version of UMFPACK code must cite the
    Copyright, this License, the Availability note, and "Used by permission."
    Permission to modify the code and to distribute modified code is granted,
    provided the Copyright, this License, and the Availability note are
    retained, and a notice that the code was modified is included.  This
    software was developed with support from the National Science Foundation,
    and is provided to you free of charge.

Availability:

    http://www.cise.ufl.edu/research/sparse/umfpack

--------------------------------------------------------------------------------

AMD Version 1.0 (Apr. 30, 2003),  Copyright (c) 2003 by Timothy A.
Davis, Patrick R. Amestoy, and Iain S. Duff.  All Rights Reserved.

AMD License:

    Your use or distribution of AMD or any modified version of
    AMD implies that you agree to this License.

    THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
    EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

    Permission is hereby granted to use or copy this program, provided
    that the Copyright, this License, and the Availability of the original
    version is retained on all copies.  User documentation of any code that
    uses AMD or any modified version of AMD code must cite the
    Copyright, this License, the Availability note, and "Used by permission."
    Permission to modify the code and to distribute modified code is granted,
    provided the Copyright, this License, and the Availability note are
    retained, and a notice that the code was modified is included.  This
    software was developed with support from the National Science Foundation,
    and is provided to you free of charge.

Availability:

    http://www.cise.ufl.edu/research/sparse/amd
*/
#include <stdio.h>
#include <string.h>
#include "Python.h"
#define PY_ARRAY_UNIQUE_SYMBOL umfpack

#include "numpy/arrayobject.h"
#include "numpy/noprefix.h"

#include "umfpack.h"
#include "spmatrix.h"
#include "ll_mat.h"

#define SYMMETRIC 1    /* Symmetric SpMatrix */
#define GENERAL   0    /* General   SpMatrix */

/***********************************************************************
 * UMFPackObject definition
 */

typedef struct UMFPackObject {
  PyObject_VAR_HEAD
  int     n;
  int     nnz;
  double *val;                  /* pointer to array of values */
  int    *row;          /* pointer to array of indices */
  int    *ind;          /* pointer to array of indices */
  void   *Numeric;
  double *Control;
  double *Info;
} UMFPackObject;

/***********************************************************************
 * UMFPackObject methods
 */

static char getlists_doc[] = "";

static PyObject *
UMFPack_getlists(UMFPackObject *self, PyObject *args) {
  PyObject *ind=NULL, *row=NULL, *val=NULL, *tupleret=NULL;
  int i;

  if (!PyArg_ParseTuple(args, ""))
    return NULL;

  ind = PyList_New(self->n+1);
  if (ind == NULL)
    goto failmemory;

  PyList_SetItem(ind, 0, Py_BuildValue("i", 0));
  for (i=1; i < (self->n + 1); i++) {
    PyList_SetItem(ind, i, Py_BuildValue("i", self->ind[i] ));
  }

  row = PyList_New(self->nnz);
  if (row == NULL)
    goto failmemory;

  for (i=0; i < self->nnz; i++) {
    PyList_SetItem(row, i, Py_BuildValue("i", self->row[i]));
  }

  val = PyList_New(self->nnz);
  if (val == NULL)
    goto failmemory;

  for (i=0; i < self->nnz; i++) {
    PyList_SetItem(val, i, Py_BuildValue("d", self->val[i]));
  }

  tupleret = PyTuple_New(3);
  if (tupleret == NULL)
    goto failmemory;

  PyTuple_SetItem(tupleret, 0, ind);
  PyTuple_SetItem(tupleret, 1, row);
  PyTuple_SetItem(tupleret, 2, val);

  return tupleret;

  failmemory:
    if (ind != NULL)
      PyObject_Del(ind);
    if (row != NULL)
      PyObject_Del(row);
    if (val != NULL)
      PyObject_Del(val);
    if (tupleret != NULL)
      PyObject_Del(tupleret);
    return PyErr_NoMemory();
}

static char lunz_doc[] = "self.lunz(): Returns number of nonzeros in factors\n";

static PyObject *UMFPack_Lunz(UMFPackObject *self, PyObject *args) {

  int     nrow, ncol, lnz, unz, nz_udiag, status;
  void   *Numeric = self->Numeric;

  /* Obtain number of nonzeros in factors */
  status = umfpack_di_get_lunz(&lnz, &unz, &nrow, &ncol, &nz_udiag, Numeric);

  if( status != UMFPACK_OK ) {
      switch(status) {
        case UMFPACK_ERROR_invalid_Numeric_object:
            PyErr_SetString(PyExc_SystemError,
                            "Get_Lunz:: Invalid Numeric object");
            return NULL;
        case UMFPACK_ERROR_argument_missing:
            PyErr_SetString(PyExc_SystemError,
                            "Get_Lunz:: Invalid arguments");
            return NULL;
      }
  }

  return Py_BuildValue("iii", lnz, unz, nz_udiag);
}

static char lu_doc[] = "self.lu()\n\
    Returns L and U factors, permutation and scaling information.\n\n\
    Use: (L, U, P, Q, R, do_recip) = self.lu().\n\
    The original matrix A is factorized into\n\
         L U = P R A Q\n\
    where L is unit lower triangular,\n\
          U is upper triangular,\n\
          P and Q are permutation matrices,\n\
          R is a row-scaling diagonal matrix such that\n\
           the i-th row of A has been divided    by R[i] if do_recip = True,\n\
           the i-th row of A has been multiplied by R[i] if do_recip = False.\n\
    \n\
    L and U are returned as ll_mat sparse matrices.\n\
    P, Q and R are returned as NumPy arrays.\
";

static PyObject *UMFPack_Lu(UMFPackObject *self, PyObject *args) {

  int     nrow, ncol, do_recip, lnz, unz, nz_udiag, status, i, j;
  int    *Lp, *Lj, *Up, *Ui;
  double *Lx, *Ux, *D;
  void   *Numeric = self->Numeric;

  PyArrayObject *P,    *Q,    *R;
  npy_int       *pval, *qval;
  npy_double                  *rval;
  npy_intp       dp[1], dq[1], dr[1];

  PyObject      *Lmat,    *Umat;
  int            ldims[2], udims[2];

  int storeZeros = 0;

  /* Obtain number of nonzeros in factors */
  status = umfpack_di_get_lunz(&lnz, &unz, &nrow, &ncol, &nz_udiag, Numeric);

  if( status != UMFPACK_OK ) {
      switch(status) {
        case UMFPACK_ERROR_invalid_Numeric_object:
            PyErr_SetString(PyExc_SystemError,
                            "Get_Lunz:: Invalid Numeric object");
            return NULL;
        case UMFPACK_ERROR_argument_missing:
            PyErr_SetString(PyExc_SystemError,
                            "Get_Lunz:: Invalid arguments");
            return NULL;
      }
  }

  /* Allocate memory to hold factors */
  Lp = (int *)calloc(nrow+1, sizeof(int));
  Lj = (int *)calloc(lnz, sizeof(int));
  Lx = (double *)calloc(lnz, sizeof(double));

  Up = (int *)calloc(ncol+1, sizeof(int));
  Ui = (int *)calloc(unz, sizeof(int));
  Ux = (double *)calloc(unz, sizeof(double));

  dp[0] = nrow; dq[0] = ncol; dr[0] = nrow;
  P = (PyArrayObject *)PyArray_SimpleNew(1, dp, NPY_INT);
  pval = (npy_int *)P->data;

  Q = (PyArrayObject *)PyArray_SimpleNew(1, dq, NPY_INT);
  qval = (npy_int *)Q->data;

  D = NULL; /* Force diagonal terms into U */

  R = (PyArrayObject *)PyArray_SimpleNew(1, dr, NPY_DOUBLE);
  rval = (npy_double *)R->data;

  /* Retrieve factors, permutation vectors and row scaling */
  status = umfpack_di_get_numeric(Lp, Lj, Lx, Up, Ui, Ux, pval, qval, D,
                                  &do_recip, rval, Numeric);

  if( status != UMFPACK_OK ) {
      switch(status) {
          case UMFPACK_ERROR_out_of_memory:
            PyErr_SetString(PyExc_SystemError, "Get_Numeric:: out of memory");
            return NULL;
          case UMFPACK_ERROR_invalid_Numeric_object:
            PyErr_SetString(PyExc_SystemError,
                            "Get_Numeric:: Numeric object is invalid");
            return NULL;
      }
  }

  /* Create LL_Mat structure for factor L and populate it */
  ldims[0] = nrow;
  ldims[1] = nrow < ncol ? nrow : ncol;
  Lmat = SpMatrix_NewLLMatObject( ldims, GENERAL, lnz, storeZeros);

  for( i=0; i<nrow+1; i++ )
    for( j=Lp[i]; j<Lp[i+1]; j++ )
        SpMatrix_LLMatSetItem((LLMatObject *)Lmat, i, Lj[j], Lx[j]);

  free( Lp ); free( Lj ); free( Lx );

  /* Create LL_Mat structure for factor U and populate it */
  udims[0] = nrow < ncol ? nrow : ncol;
  udims[1] = ncol;
  Umat = SpMatrix_NewLLMatObject( udims, GENERAL, unz, storeZeros);

  for( j=0; j<ncol+1; j++ )
      for( i=Up[j]; i<Up[j+1]; i++ )
          SpMatrix_LLMatSetItem((LLMatObject *)Umat, Ui[i], j, Ux[i]);

  free( Up ); free( Ui ); free( Ux );

  /* Return output data */
  return Py_BuildValue( "OOOOOi", Lmat, Umat, P, Q, R, do_recip );
}

static char solve_doc[] = "self.solve(b, x, method, irsteps)\n\
\n\
solves linear system of equations.\n\
\n\
parameters\n\
----------\n\
\n\
b        array, right hand side of equation\n\
x        array, solution vector\n\n\
method   'UMFPACK_A': solve A   * x == b\n\
         'UMFPACK_At': solve A^T * x == b\n\
         'UMFPACK_Pt_L': solve P^T * L * x == b\n\
         'UMFPACK_L': solve L * x == b\n\
         'UMFPACK_Lt_P': solve L^T * P * x == b\n\
         'UMFPACK_Lt': solve L^T * x == b\n\
         'UMFPACK_U_Qt': solve U * Q^T * x == b\n\
         'UMFPACK_U': solve U * x == b\n\
         'UMFPACK_Q_Ut': solve Q * U^T * x == b\n\
         'UMFPACK_Ut': solve U^T * x == b\n\
\n\
irsteps  Number of iterative refinement steps to attempt.\n\
";

static PyObject *
UMFPack_solve(UMFPackObject *self, PyObject *args) {
  PyArrayObject *b, *x;
  char *method = "UMFPACK_A";
  int res, sys, irsteps = -1, status;

  res = PyArg_ParseTuple(args, "O!O!|si",
                               &PyArray_Type, &b,
                               &PyArray_Type, &x,
                               &method, &irsteps);
  if( !res )
      return NULL;

  SPMATRIX_CHECK_ARR_DIM_SIZE(b, 1, self->n);
  SPMATRIX_CHECK_ARR_DIM_SIZE(x, 1, self->n);

  if (strcmp(method, "UMFPACK_A") == 0)
      sys = UMFPACK_A;
  else if (strcmp(method, "UMFPACK_At") == 0)
      sys = UMFPACK_At;
  else if (strcmp(method, "UMFPACK_Pt_L") == 0)
      sys = UMFPACK_Pt_L;
  else if (strcmp(method, "UMFPACK_L") == 0)
      sys = UMFPACK_L;
  else if (strcmp(method, "UMFPACK_Lt_P") == 0)
      sys = UMFPACK_Lt_P;
  else if (strcmp(method, "UMFPACK_Lt") == 0)
      sys = UMFPACK_Lt;
  else if (strcmp(method, "UMFPACK_U_Qt") == 0)
      sys = UMFPACK_U_Qt;
  else if (strcmp(method, "UMFPACK_U") == 0)
      sys = UMFPACK_U;
  else if (strcmp(method, "UMFPACK_Q_Ut") == 0)
      sys = UMFPACK_Q_Ut;
  else if (strcmp(method, "UMFPACK_Ut") == 0)
      sys = UMFPACK_Ut;
  else {
      PyErr_SetString(PyExc_ValueError, "method");
      return NULL;
  }

  if( irsteps != -1 )
      self->Control[UMFPACK_IRSTEP] = irsteps;

  status = umfpack_di_solve(sys, self->ind, self->row, self->val,
                            (double *)x->data, (double *)b->data,
                            self->Numeric, self->Control, self->Info);

  switch( status ) {
  case UMFPACK_WARNING_singular_matrix:
    PyErr_SetString(PyExc_SystemError, "singular matrix");
    return NULL;
  case UMFPACK_ERROR_out_of_memory:
    PyErr_SetString(PyExc_SystemError, "insufficient memory");
    return NULL;
  case UMFPACK_ERROR_argument_missing:
    PyErr_SetString(PyExc_SystemError, "one or more argument missing");
    return NULL;
  case UMFPACK_ERROR_invalid_system:
    PyErr_SetString(PyExc_SystemError, "matrix is not square");
    return NULL;
  case UMFPACK_ERROR_invalid_Numeric_object:
    PyErr_SetString(PyExc_SystemError, "invalid Numeric object");
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

/** table of object methods
 */
PyMethodDef UMFPack_methods[] = {
  {"solve", (PyCFunction)UMFPack_solve, METH_VARARGS, solve_doc},
  {"getlists", (PyCFunction)UMFPack_getlists, METH_VARARGS, getlists_doc},
  {"lunz", (PyCFunction)UMFPack_Lunz, METH_VARARGS, lunz_doc},
  {"lu", (PyCFunction)UMFPack_Lu, METH_VARARGS, lu_doc},
  {NULL, NULL}          /* sentinel */
};


/***********************************************************************
 * UMFPackType methods
 */

static void
UMFPack_dealloc(UMFPackObject *self)
{
  umfpack_di_free_numeric(&(self->Numeric));
  PyMem_Free(self->val);
  PyMem_Free(self->row);
  PyMem_Free(self->ind);
  PyMem_Free(self->Control);
  PyMem_Free(self->Info);
  PyObject_Del(self);
}

static PyObject *
UMFPack_getattr(UMFPackObject *self, char *name)
{
  if (strcmp(name, "shape") == 0)
    return Py_BuildValue("(i,i)", self->n, self->n);
  if (strcmp(name, "nnz") == 0)
    return Py_BuildValue("i", self->nnz);
  if (strcmp(name, "__members__") == 0) {
    char *members[] = {"shape", "nnz"};
    int i;

    PyObject *list = PyList_New(sizeof(members)/sizeof(char *));
    if (list != NULL) {
      for (i = 0; i < sizeof(members)/sizeof(char *); i ++)
    PyList_SetItem(list, i, PyString_FromString(members[i]));
      if (PyErr_Occurred()) {
    Py_DECREF(list);
    list = NULL;
      }
    }
    return list;
  }
  return Py_FindMethod(UMFPack_methods, (PyObject *)self, name);
}


/***********************************************************************
 * UMFPackType structure
 */

PyTypeObject UMFPackType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "umfpack_context",
  sizeof(UMFPackObject),
  0,
  (destructor)UMFPack_dealloc,   /* tp_dealloc */
  0,                /* tp_print */
  (getattrfunc)UMFPack_getattr,  /* tp_getattr */
  0,                /* tp_setattr */
  0,                /* tp_compare */
  0,                /* tp_repr */
  0,                /* tp_as_number*/
  0,                /* tp_as_sequence*/
  0,                /* tp_as_mapping*/
  0,                /* tp_hash */
};

/***********************************************************************
 * Object construction functions
 */

/***********************************************************************
 * Module functions
 */

static void sortcol(int *row, double *val, int n)
{
  int i, j, tmpint;
  double tmpdbl;

  for (i=n-1; i > 0; i--) {

    for (j=0; j < i; j++) {
      if (row[j] > row[j+1]) {

        /* swap row indices */
        tmpint = row[j+1];
        row[j+1] = row[j];
        row[j] = tmpint;

        /* swap values */
        tmpdbl = val[j+1];
        val[j+1] = val[j];
        val[j] = tmpdbl;

      }
    }
  }

  return;
}

// newUMFPackObject(LLMatObject *matrix, int strategy, double tol2by2, int scale,
//                   double tolpivot, double tolsympivot)
static PyObject *
newUMFPackObject(LLMatObject *matrix, int strategy, int scale,
                  double tolpivot, double tolsympivot)
{
  UMFPackObject *self;
  void *Symbolic;
  struct llColIndex *colidx;
  int i, validx, curridx, status, exitCode=0;

  if (SpMatrix_LLMatBuildColIndex(&colidx, matrix, 1) == 1)
    return NULL;

  self = PyObject_New(UMFPackObject, &UMFPackType);
  if (self == NULL)
    return PyErr_NoMemory();

  /* initialize pointers to arrays */
  self->val = NULL;
  self->row = NULL;
  self->ind = NULL;
  self->Control = NULL;
  self->Info = NULL;
  self->Numeric = NULL;

  self->n = matrix->dim[0];
  if (matrix->issym)
    self->nnz = 2 * colidx->nzLo + colidx->nzDiag + 2 * colidx->nzUp;
  else
    self->nnz = matrix->nnz;

  self->val = (double *)PyMem_Malloc(self->nnz * sizeof(double));
  if (self->val == NULL) {
    exitCode = -1;
    goto freememory;
  }

  self->row = (int *)PyMem_Malloc(self->nnz * sizeof(int));
  if (self->row == NULL) {
    exitCode = -1;
    goto freememory;
  }
  self->ind = (int *)PyMem_Malloc((self->n + 1) * sizeof(int) );
  if (self->ind == NULL) {
    exitCode = -1;
    goto freememory;
  }
  self->ind[self->n] = self->nnz;

  self->Control = (double *)PyMem_Malloc(UMFPACK_CONTROL * sizeof(double));
  if (self->Control == NULL) {
    exitCode = -1;
    goto freememory;
  }

  umfpack_di_defaults(self->Control);

  if (strategy != -1)
    self->Control[UMFPACK_STRATEGY] = strategy;

  //if (tol2by2 != -1)
  //  self->Control[UMFPACK_2BY2_TOLERANCE] = tol2by2;

  if (scale != -1)
    self->Control[UMFPACK_SCALE] = scale;

  if (tolpivot != -1)
    self->Control[UMFPACK_PIVOT_TOLERANCE] = tolpivot;

  if (tolsympivot != -1)
    self->Control[UMFPACK_SYM_PIVOT_TOLERANCE] = tolsympivot;

  self->Info = (double *)PyMem_Malloc(UMFPACK_INFO * sizeof(double));
  if (self->Info == NULL) {
    exitCode = -1;
    goto freememory;
  }

  validx = 0;
  for (i=0; i < self->n; i++) {
    self->ind[i] = validx;
    for (curridx=colidx->root[i]; curridx != -1; curridx=colidx->link[curridx]){
      self->val[validx] = matrix->val[curridx];
      self->row[validx++] = colidx->row[curridx];
    }

    if (matrix->issym) {
      for (curridx=matrix->root[i]; curridx != -1; curridx=matrix->link[curridx]){
        if (i != matrix->col[curridx]) {
          self->val[validx] = matrix->val[curridx];
          self->row[validx++] = matrix->col[curridx];
        }
      }
    }

    sortcol(&(self->row[self->ind[i]]), &(self->val[self->ind[i]]), validx - self->ind[i]);
  }

  status = umfpack_di_symbolic(self->n, self->n, self->ind, self->row, self->val, &Symbolic, self->Control, self->Info);

  switch(status) {
    case UMFPACK_ERROR_n_nonpositive:
      PyErr_SetString(PyExc_SystemError, "n must be positive");
      exitCode = -2;
      goto freememory;
    case UMFPACK_ERROR_out_of_memory:
      PyErr_SetString(PyExc_SystemError, "insufficient memory");
      exitCode = -2;
      goto freememory;
    case UMFPACK_ERROR_argument_missing:
      PyErr_SetString(PyExc_SystemError, "one or more argument missing");
      exitCode = -2;
      goto freememory;
    case UMFPACK_ERROR_invalid_matrix:
      PyErr_SetString(PyExc_SystemError, "invalid matrix");
      exitCode = -2;
      goto freememory;
    case UMFPACK_ERROR_internal_error:
      PyErr_SetString(PyExc_SystemError, "bug in umfpack");
      exitCode = -2;
      goto freememory;
  }

  status = umfpack_di_numeric(self->ind, self->row, self->val, Symbolic, &(self->Numeric), self->Control, self->Info);

  switch(status) {

    case UMFPACK_WARNING_singular_matrix:
      PyErr_SetString(PyExc_SystemError, "singular matrix");
      exitCode = -2;
      goto freememory;
    case UMFPACK_ERROR_out_of_memory:
      PyErr_SetString(PyExc_SystemError, "insufficient memory");
      exitCode = -2;
      goto freememory;
    case UMFPACK_ERROR_argument_missing:
      PyErr_SetString(PyExc_SystemError, "one or more argument missing");
      exitCode = -2;
      goto freememory;
    case UMFPACK_ERROR_invalid_Symbolic_object:
      PyErr_SetString(PyExc_SystemError, "invalid symbolic object");
      exitCode = -2;
      goto freememory;
    case UMFPACK_ERROR_different_pattern:
      PyErr_SetString(PyExc_SystemError, "matrix has changed since last factorization");
      exitCode = -2;
      goto freememory;
  }

  umfpack_di_free_symbolic(&Symbolic);

  return (PyObject *)self;

  freememory:
    SpMatrix_LLMatDestroyColIndex(&colidx);
    if (self->val != NULL)
      PyMem_Free(self->val);
    if (self->row != NULL)
      PyMem_Free(self->row);
    if (self->ind != NULL)
      PyMem_Free(self->ind);
    if (self->Control != NULL)
      PyMem_Free(self->Control);
    if (self->Info != NULL)
      PyMem_Free(self->Info);
    PyObject_Del(self);

    if (exitCode == -1)
        return PyErr_NoMemory();
    else
        return NULL;

}

static char factorize_doc[] = "factorize(A, ...)\n\
\n\
performs a factorization of the sparse matrix A (which is a ll_mat object) and \n\
return a umfpack_context object.\n\
\n\
arguments\n\
---------\n\
\n\
A    spmatrix.ll_mat object.\n\
     Matrix to be factorized\n\
\n\
additional keyword arguments:\n\
-----------------------------\n\
\n\
strategy            determines what kind of ordering and\n\
                    pivoting strategy that UMFPACK should use.\n\
                    \"UMFPACK_STRATEGY_AUTO\"\n\
                    \"UMFPACK_STRATEGY_UNSYMMETRIC\"\n\
                    \"UMFPACK_STRATEGY_SYMMETRIC\"\n\
                    \"UMFPACK_STRATEGY_2BY2\"\n\
\n\
tol2by2             tolerance for the 2 by 2 strategy\n\
\n\
scale               scaling UMFPACK would use\n\
                    \"UMFPACK_SCALE_NONE\"\n\
                    \"UMFPACK_SCALE_SUM\"\n\
                    \"UMFPACK_SCALE_MAX\"\n\
\n\
tolpivot            relative pivot tolerance for threshold partial\n\
                    pivoting with row interchanges.\n\
\n\
tolsympivot         if diagonal pivoting is attempted then this\n\
                    parameter is used to control when the diagonal\n\
                    is selected in a given pivot column.\n\
";

static PyObject *
factorize(PyObject *self, PyObject *args, PyObject *keywds) {
  LLMatObject *matrix;
  char *strategy="UMFPACK_STRATEGY_AUTO";
  //double tol2by2 = 0.1;
  char *scale="UMFPACK_SCALE_SUM";
  double tolpivot = 0.1;
  double tolsympivot = 0.0;
  int res;
  int strategyval = UMFPACK_STRATEGY_AUTO, scaleval = UMFPACK_SCALE_SUM;

  //static char *kwlist[] = {"", "strategy", "tol2by2", "scale", "tolpivot", "tolsympivot", NULL};
  static char *kwlist[] = {"", "strategy", "scale", "tolpivot", "tolsympivot", NULL};

  /*
   res = PyArg_ParseTupleAndKeywords(args, keywds, "O!|sdsdd", kwlist,
                    &LLMatType, &matrix,
                    &strategy,
                    &tol2by2,
                    &scale,
                    &tolpivot,
                    &tolsympivot);
                    */
  res = PyArg_ParseTupleAndKeywords(args, keywds, "O!|ssdd", kwlist,
                    &LLMatType, &matrix,
                    &strategy,
                    &scale,
                    &tolpivot,
                    &tolsympivot);
  if (!res)
    return NULL;

  if (strcmp("UMFPACK_STRATEGY_AUTO", strategy) == 0)
    strategyval = UMFPACK_STRATEGY_AUTO;
  else if (strcmp("UMFPACK_STRATEGY_UNSYMMETRIC", strategy) == 0)
    strategyval = UMFPACK_STRATEGY_UNSYMMETRIC;
  else if (strcmp("UMFPACK_STRATEGY_SYMMETRIC", strategy) == 0)
    strategyval = UMFPACK_STRATEGY_SYMMETRIC;
  // else if (strcmp("UMFPACK_STRATEGY_2BY2", strategy) == 0)
  //   strategyval = UMFPACK_STRATEGY_2BY2;

  if (strcmp("UMFPACK_SCALE_NONE", scale) == 0)
    scaleval = UMFPACK_SCALE_NONE;
  else if (strcmp("UMFPACK_SCALE_SUM", scale) == 0)
    scaleval = UMFPACK_SCALE_SUM;
  if (strcmp("UMFPACK_SCALE_MAX", scale) == 0)
    scaleval = UMFPACK_SCALE_MAX;

  // return newUMFPackObject(matrix, strategyval, tol2by2, scaleval, tolpivot, tolsympivot);
  return newUMFPackObject(matrix, strategyval, scaleval, tolpivot, tolsympivot);
}


/** table of module functions
 */
static PyMethodDef precon_methods[] = {
  {"factorize", (PyCFunction)factorize, METH_VARARGS|METH_KEYWORDS, factorize_doc},
  {NULL, NULL}  /* sentinel */
};


DL_EXPORT(void)
initumfpack(void)
{
  PyObject *m, *d;

  UMFPackType.ob_type = &PyType_Type;

  m = Py_InitModule("umfpack", precon_methods);
  d = PyModule_GetDict(m);

  PyDict_SetItemString(d, "UMFPackType", (PyObject *)&UMFPackType);

  /* initialize Numeric array module */
  import_array();
  /* initialize spmatrix module */
  import_spmatrix();

  /* No need to check the error here, the caller will do that */
}
