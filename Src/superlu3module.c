#include <stdio.h>
#include "Python.h"
#define PY_ARRAY_UNIQUE_SYMBOL superlu

#if NUMPY
#include "numpy/arrayobject.h"
#include "numpy/noprefix.h"
#else
#include "Numeric/arrayobject.h"
#endif

#include "slu_ddefs.h"
#include "slu_util.h"

#include "pysparse/spmatrix.h"

/*********************************************************************** 
 * static variable
 */


/*********************************************************************** 
 * SuperLUObject definition
 */

typedef struct SuperLUObject {
  PyObject_VAR_HEAD
  int n;
  SuperMatrix L;
  SuperMatrix U;
  int *perm_r;
  int *perm_c;
  SuperLUStat_t stat;
  int StatInit_done;  /* flag showing whether StatInit was already called */
  superlu_options_t options;
} SuperLUObject;

/*********************************************************************** 
 * SuperLUObject methods
 */

static char solve_doc[] = "self.solve(b, x, trans)\n\
\n\
solves linear system of equations with one or sereral right hand sides.\n\
\n\
parameters\n\
----------\n\
\n\
b        array, right hand side(s) of equation\n\
x        array, solution vector(s)\n\
trans    'N': solve A   * x == b\n\
         'T': solve A^T * x == b\n\
         (optional, default value 'N')\n\
";

static PyObject *SuperLU_solve(SuperLUObject *self, PyObject *args) {
  PyArrayObject *b, *x;
  SuperMatrix B;
  char trans = 'N';
  trans_t Trans;
  int i, info = 0;

  if (!PyArg_ParseTuple(args, "O!O!|c", 
			&PyArray_Type, &b, 
			&PyArray_Type, &x,
			&trans))
    return NULL;
  
  SPMATRIX_CHECK_ARR_DIM_SIZE(b, 1, self->n);
  SPMATRIX_CHECK_ARR_DIM_SIZE(x, 1, self->n);

  /* solve transposed system:
   * matrix was passed row-wise instead of column-wise */
  if (trans == 'n' || trans == 'N')
    Trans = TRANS; //trans = 'T';
  else if (trans == 't' || trans == 'T')
    Trans = NOTRANS; //trans = 'N';
  else {
    PyErr_SetString(PyExc_ValueError, "trans");
    return NULL;
  }

  /* copy b to x */
  for (i = 0; i < self->n; i ++)
    ((double *)x->data)[i] = ((double *)b->data)[i];

  /* Create data structure for right hand side */
  dCreate_Dense_Matrix(&B, self->n, 1, (double *)x->data, self->n,
                       SLU_DN, SLU_D, SLU_GE);

  /* Solve the system, overwriting vector x. */
  dgstrs(Trans, &self->L, &self->U, self->perm_c, self->perm_r, &B,
         &self->stat, &info);

  /* free memory */
  Destroy_SuperMatrix_Store(&B);

  if (info) {
    PyErr_SetString(PyExc_SystemError,
                    "dgstrs was called with invalid arguments");
    return NULL;
  } else {
    Py_INCREF(Py_None); 
    return Py_None;
  }
}

/** table of object methods
 */
PyMethodDef SuperLU_special_methods[] = {
  {"solve", (PyCFunction)SuperLU_solve, METH_VARARGS, solve_doc},
  {NULL,    NULL,                       0,            NULL     }  /* sentinel */
};


/*********************************************************************** 
 * SuperLUType methods
 */

static void SuperLU_dealloc(SuperLUObject *self) {
  SUPERLU_FREE(self->perm_r);
  SUPERLU_FREE(self->perm_c);
  Destroy_SuperNode_Matrix(&self->L);
  Destroy_CompCol_Matrix(&self->U);
  StatFree(&self->stat);
  PyObject_Del(self);
}

static PyObject *SuperLU_getattr(SuperLUObject *self, char *name) {
  if (strcmp(name, "shape") == 0)
    return Py_BuildValue("(i,i)", self->n, self->n);
  if (strcmp(name, "nnz") == 0)
    return Py_BuildValue("i", ((SCformat *)self->L.Store)->nnz + ((SCformat *)self->U.Store)->nnz);
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
  return Py_FindMethod(SuperLU_special_methods, (PyObject *)self, name);
}


/***********************************************************************
 * SuperLUType structure
 */

PyTypeObject SuperLUType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "superlu_context",
  sizeof(SuperLUObject),
  0,
  (destructor)SuperLU_dealloc,   /* tp_dealloc */
  0,				/* tp_print */
  (getattrfunc)SuperLU_getattr,  /* tp_getattr */
  0,                           /* tp_setattr */
  0,                           /* tp_compare */
  0,                           /* tp_repr */
  0,                           /* tp_as_number*/
  0,                           /* tp_as_sequence*/
  0,                           /* tp_as_mapping*/
  0,                           /* tp_hash */
  0,                           /* tp_call*/
  0,                           /* tp_str*/
  0,                           /* tp_getattro*/
  0,                           /* tp_setattro*/
  0,                           /* tp_as_buffer*/
  Py_TPFLAGS_DEFAULT,          /* tp_flags*/
  "SuperLU Context Object",    /* tp_doc */
};

/*********************************************************************** 
 * Object construction functions
 */

/*********************************************************************** 
 * Module functions
 */

static PyObject *newSuperLUObject(int n, CSRMatObject *matrix, PyObject *issym,
                                  double diag_pivot_thresh, double drop_tol,
                                  int relax, int panel_size, int permc_spec) {
  SuperLUObject *self;
  SuperMatrix A;      /* A in NC format used by the factorization routine. */
  SuperMatrix AC;     /* Matrix postmultiplied by Pc */
  mem_usage_t mem_usage;
  int *etree;
  int info = 0;

  //printf("Creating new object\n");
  
  /* Create SuperLUObject */
  self = PyObject_New(SuperLUObject, &SuperLUType);
  if (self == NULL)
    return NULL; // PyErr_NoMemory();
  self->n = n;
  self->perm_r = NULL;
  self->perm_c = NULL;
  self->StatInit_done = 0;

  //printf("Setting options (permc_spec = %d)\n", permc_spec);

  set_default_options(&self->options);
  switch( permc_spec ) {
  case 0:
    self->options.ColPerm = NATURAL;
    break;
  case 1:
    self->options.ColPerm = MMD_ATA;
    break;
  case 2:
    self->options.ColPerm = MMD_AT_PLUS_A;
    break;
  case 3:
    self->options.ColPerm = COLAMD;
    break;
  default:
    //self->options.ColPerm = MY_PERMC;
    break;
  }
  self->options.DiagPivotThresh = diag_pivot_thresh;
  self->options.PrintStat = YES;
  self->options.SymmetricMode = (issym == Py_True) ? YES : NO;

  //printf("Setting up stats\n");

  /* Make sure StatInit is only called once */
  if(! self->StatInit_done) {
    StatInit(&self->stat);    // StatInit(panel_size, relax);
    self->StatInit_done = 1;
  }

  //printf("Creating matrix structure\n");
  /* Create matrix structure ; indicate that A is stored by rows (SLU_NR) */
  dCreate_CompCol_Matrix(&A, n, n, matrix->nnz,
                         matrix->val, matrix->col, matrix->ind,
                         SLU_NR, SLU_D, (issym == Py_True) ? SLU_SYL : SLU_GE);

  //dPrint_CompCol_Matrix("A", &A);

  //printf("Allocating room for permutations\n");
  etree = intMalloc(n);
  self->perm_r = intMalloc(n);
  self->perm_c = intMalloc(n);
  if (etree == NULL || self->perm_r == NULL || self->perm_c == NULL) {
    PyErr_NoMemory();
    goto fail;
  }

  //printf("Obtaining permutation\n");
  /* Obtain and apply column permutation */
  get_perm_c(permc_spec, &A, self->perm_c);
  sp_preorder(&self->options, &A, self->perm_c, etree, &AC);

  //print_int_vec("\nperm_c", n, self->perm_c);

  //printf("Factorizing\n");
  /* Perform factorization (perm_c and perm_r are swapped because our matrix
   * is stored in compressed-row format and not in compressed-column format! */
  dgstrf(&self->options, &AC, drop_tol, relax, panel_size,
         etree, NULL, 0, self->perm_c, self->perm_r,
         &self->L, &self->U, &self->stat, &info);

  //dPrint_SuperNode_Matrix("L", &self->L);
  //dPrint_CompCol_Matrix("U", &self->U);

  //printf("Freeing memory\n");
  /* free memory */
  SUPERLU_FREE(etree);
  Destroy_SuperMatrix_Store(&A); /* arrays of input matrix will not be freed */
  Destroy_CompCol_Permuted(&AC);

  if (info) {
    if (info < 0)
      PyErr_SetString(PyExc_SystemError,
                      "dgstrf was called with invalid arguments");
    else {
      if (info <= n) 
	PyErr_SetString(PyExc_RuntimeError, "Factor is exactly singular");
      else
	PyErr_NoMemory();
    }
    goto fail;
  }
  
  return (PyObject *)self;

 fail:
  PyMem_Del(self->perm_r);
  PyMem_Del(self->perm_c);
  StatFree(&self->stat);
  PyObject_Del(self);
  return NULL;
}

static char factorize_doc[] = "factorize(A, ...)\n\
\n\
performs a factorization of the sparse matrix A (which is a csr_mat object) and \n\
return a superlu_context object.\n\
\n\
arguments\n\
---------\n\
\n\
A    spmatrix.csr_mat object.\n\
     Matrix to be factorized\n\
\n\
additional keyword arguments:\n\
-----------------------------\n\
\n\
diag_pivot_thresh   threshhold for partial pivoting.\n\
                    0.0 <= diag_pivot_thresh <= 1.0\n\
                    0.0 corresponds to no pivoting\n\
                    1.0 corresponds to partial pivoting\n\
                    (default: 1.0)\n\
\n\
drop_tol            drop tolerance parameter\n\
                    0.0 <= drop_tol <= 1.0\n\
                    0.0 corresponds to exact factorization\n\
                    CAUTION: the drop_tol is not implemented in SuperLU 2.0\n\
                    (default: 0.0)\n\
\n\
relax               to control degree of relaxing supernodes\n\
                    (default: 1)\n\
\n\
panel_size          a panel consist of at most panel_size consecutive columns.\n\
                    (default: 10)\n\
\n\
permc_spec          specifies the matrix ordering used for the factorization\n\
                    0: natural ordering\n\
                    1: MMD applied to the structure of A^T * A\n\
                    2: MMD applied to the structure of A^T + A\n\
                    3: COLAMD, approximate minimum degree column ordering\n\
                    (default: 2)\n\
";

/* ========================================================================== */

static PyObject *factorize(PyObject *self, PyObject *args, PyObject *keywds) {
  int n;			/* dimension of matrix */

  /* default value for SuperLU parameters*/
  double diag_pivot_thresh = 1.0;
  double drop_tol = 0.0;
  PyObject *symmetric = Py_False;
  int relax = 1;
  int panel_size = 10;
  int permc_spec = 2;
  
  PyObject *matrix;
  static char *kwlist[] = {"", "symmetric", "diag_pivot_thresh", "drop_tol",
                           "relax", "panel_size", "permc_spec", NULL};

  int res = PyArg_ParseTupleAndKeywords(args, keywds, "O|Oddiii", kwlist, 
					&matrix,
                                        &symmetric,
					&diag_pivot_thresh,
					&drop_tol,
					&relax,
					&panel_size,
					&permc_spec);
  if (!res) return NULL;

  /* check shape of matrix object */
  if (SpMatrix_GetOrder(matrix, &n)) return NULL;

  return newSuperLUObject(n, (CSRMatObject *)matrix, symmetric,
                          diag_pivot_thresh, drop_tol, relax, panel_size,
                          permc_spec);
}

/* ========================================================================== */

/*
 *          D e f i n i t i o n   o f   S u p e r L U   m e t h o d s
 */

/* ========================================================================== */

static PyMethodDef SuperLU_methods[] = {
  {"factorize", (PyCFunction)factorize, METH_VARARGS|METH_KEYWORDS, factorize_doc},
  {NULL,        NULL,      0,                          NULL         }// sentinel
};

/* ========================================================================== */

DL_EXPORT(void) initsuperlu(void) {
  PyObject *m, *d;
  
  //SuperLUType.ob_type = &PyType_Type;
  if( PyType_Ready( &SuperLUType ) < 0 ) return;

  m = Py_InitModule3("superlu", SuperLU_methods, "Python interface to SuperLU");
  d = PyModule_GetDict(m);
  PyDict_SetItemString(d, "SuperLUType", (PyObject *)&SuperLUType);

  import_array();     /* Initialize the NumPy module */
  import_spmatrix();  /* Initialize PySparse module */

  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("Unable to initialize module superlu");

  return;
}
