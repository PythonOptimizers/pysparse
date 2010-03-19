#include "Python.h"
#include <math.h>
#include "pysparse/mmio.h"
#include "pysparse/ll_mat.h"
#include "pysparse/csr_mat.h"
#include "pysparse/sss_mat.h"

#define SPMATRIX_MODULE
//#include "pysparse/spmatrix.h"

#if NUMPY
    #include "numpy/arrayobject.h"
    #include "numpy/noprefix.h"
#else
    #include "Numeric/arrayobject.h"
#endif

#define PY_ARRAY_UNIQUE_SYMBOL spmatrix
//#include "Numeric/arrayobject.h"

#define INCREASE_FACTOR   1.5 // increase rate for reallocation of ll_mat arrays
#define PPRINT_ROW_THRESH 500 // row threshold for choosing print format
#define PPRINT_COL_THRESH 20  // column threshold for choosing print format
#define OPT_MATMATMUL     1   // enable optimised matrix-matrix-multiplication

#define PYSP_MAX(A,B) ( (A) > (B) ? (A) : (B) )

// Iterator object
#define SLICE 0
#define LIST  1
#define ARRAY 2
typedef struct PysparseIterator {

  int type;                        // Slice, list or array
  PyObject *object;                // Python slice, list or Numpy array
  int (*init)(void *);             // Iterator initialization function
  int (*notDone)(void *);          // Determines whether iterations are finished
  void (*next)(void *);            // Move to next iteration
  PyObject* (*data)(void *);       // Get current iteration data
  long (*size)(void *);            // Get total number of iterations
  long counter;
  long length;
  long start;
  long step;
  long stop;

} PysparseIterator;

// Prototypes

static PyObject *LLMat_Find(LLMatObject *self, PyObject *args);
PysparseIterator* PysparseIterator_Create( int type, PyObject *object );
void PysparseIterator_Destroy( PysparseIterator **iterator );
int PysparseIterator_List_Init( void *self );
int PysparseIterator_List_NotDone( void *self );
void PysparseIterator_List_Next( void *self );
PyObject *PysparseIterator_List_Data( void *self );
long PysparseIterator_List_Size( void *self );

// Iterator creation function: No argument check for now
PysparseIterator* PysparseIterator_Create( int type, PyObject *object ) {
  PysparseIterator *iterator;
  iterator = calloc( 1, sizeof(PysparseIterator) );
  if( iterator == NULL ) return NULL;
  iterator->type = type;
  iterator->object = object;
  iterator->init = PysparseIterator_List_Init;
  iterator->notDone = PysparseIterator_List_NotDone;
  iterator->next = PysparseIterator_List_Next;
  iterator->data = PysparseIterator_List_Data;
  iterator->size = PysparseIterator_List_Size;
  return iterator;
}

// Iterator destruction function: Let garbage collector take care of PyObjects
// Call as: PysparseIterator_Destroy( &iterator );
void PysparseIterator_Destroy( PysparseIterator **iterator ) {
  if( *iterator ) {
    free(*iterator);
    *iterator = NULL;
  }
  return;
}

// Initialization of an iterator around a Python list
int PysparseIterator_List_Init( void *self ) {

  PysparseIterator *iterator = (PysparseIterator *)self;
  if( !PyList_Check( iterator->object ) ) return -1;
  iterator->counter = 0;
  iterator->length = (long)PyList_Size( (PyObject *)(iterator->object) );
  iterator->start = 0;
  iterator->step = 1;
  iterator->stop = iterator->length;
  return 0;
}

// Define iterator functions for a Python list
int PysparseIterator_List_NotDone( void *self ) {
  PysparseIterator *iterator = (PysparseIterator *)self;
  return (iterator->counter < iterator->length);
}

void PysparseIterator_List_Next( void *self ) {
  PysparseIterator *iterator = (PysparseIterator *)self;
  iterator->counter++;
  return;
}

PyObject *PysparseIterator_List_Data( void *self ) {
  PysparseIterator *iterator = (PysparseIterator *)self;
  return (PyObject *)PyList_GetItem(iterator->object,
                                    (Py_ssize_t)(iterator->counter));
}

long PysparseIterator_List_Size( void *self ) {
  PysparseIterator *iterator = (PysparseIterator *)self;
  return iterator->length;
}

// Define iterator functions for a Numpy array




/******************************************************************************/
/*  Routines for building data structure for  columnwise traversal of ll_mat  */
/*  sparse matrix objects.                                                    */
/******************************************************************************/

/*
 *  SpMatrix_LLMatBuildColIndex: build data structure for columnwise traversal
 *
 *    build a (root,row,link) linked-list data structure which links the
 *    entries of each column of self.
 *
 *    if includeDiagonal is zero the diagonal elements of self are
 *    not included in the linked-list data structure.
 */

static int 
SpMatrix_LLMatBuildColIndex(struct llColIndex **idx, LLMatObject *self,
                            int includeDiagonal) {

  int i, j, k;

  if(!(*idx = (struct llColIndex*)malloc(sizeof(struct llColIndex)))) goto fail;
  
  /* Allocate (link,row,root) arrays. Can we do better than nalloc ??? */
  if( !( (*idx)->link = PyMem_New(int, self->nalloc) ) ) goto fail;
  if( !( (*idx)->row  = PyMem_New(int, self->nalloc) ) ) goto fail;
  if( !( (*idx)->root = PyMem_New(int, self->dim[1]) ) ) goto fail;

  /* Initialize root arrays */
  for (i = 0; i < self->dim[1]; i ++) (*idx)->root[i] = -1;

  /* Scan matrix from bottom up so the resulting lists
   * are sorted by ascending row */
  (*idx)->nzLo = 0; (*idx)->nzDiag = 0; (*idx)->nzUp = 0;
  for (i = self->dim[0] - 1; i >= 0; i--) {
    k = self->root[i];
    while (k != -1) {
      j = self->col[k];
      if (i > j)
	(*idx)->nzLo++;
      else if (i == j)
	(*idx)->nzDiag++;
      else
	(*idx)->nzUp++;
      if (includeDiagonal || i != j) {
	(*idx)->link[k] = (*idx)->root[j];
	(*idx)->root[j] = k;
	(*idx)->row[k] = i;
      }
      k = self->link[k];
    }
  }
  return 0;

 fail:
  if (*idx != NULL) {
    PyMem_Del((*idx)->link);    
    PyMem_Del((*idx)->row);    
    PyMem_Del((*idx)->root);
    free(*idx);
    *idx = NULL;
  }
  PyErr_NoMemory();
  return 1;
}

/*
 *  Free memory of (root,row,link) linked-list data
 */

static void
SpMatrix_LLMatDestroyColIndex(struct llColIndex **idx) {

  if (*idx != NULL) {
    PyMem_Del((*idx)->link);    
    PyMem_Del((*idx)->row); 
    PyMem_Del((*idx)->root);
    free(*idx);
    *idx = NULL;
  }
}  

/******************************************************************************/
/*  Routines for setting, updating and reading entries of ll_mat objects      */
/******************************************************************************/

/*
 *  Return matrix entry a[i,j] as a double value
 */

static double
SpMatrix_LLMatGetItem(LLMatObject *a, int i, int j) {

  int k, t;
  
  /*
  while( i < 0 )
    i += a->dim[0];

  while( j < 0 )
    j += a->dim[1];

  if( i >= a->dim[0] || j >= a->dim[1] ) {
    PyErr_SetString(PyExc_IndexError, "Index out of bounds");
    return 0.0;
  }
pwd  */

  if(i < 0 || i >= a->dim[0] || j < 0 || j >= a->dim[1]) {
    PyErr_SetString(PyExc_IndexError, "indices out of range");
    return 0.0;
  }

  if (a->issym && i < j) {
    t = i; i = j; j = t;
  }

  k = a->root[i];
  while (k != -1) {
    if (a->col[k] == j)
      return a->val[k];
    k = a->link[k];
  }
  return 0.0;
}

/*
 *  Set matrix entry: a[i,j] = x
 */

static int 
SpMatrix_LLMatSetItem(LLMatObject *a, int i, int j, double x) {

  void *temp;
  int k, new_elem, last, col;

  if (a->issym && i < j) {
    PyErr_SetString(PyExc_IndexError, //SpMatrix_ErrorObject, 
		    "write operation to upper triangle of symmetric matrix");
    return -1;
  } 

  //printf("LLMatSetItem: matrix dimensions=(%d,%d)\n", a->dim[0], a->dim[1]);
  //printf("LLMatSetItem: receiving row=%d, col=%d, val=%g\n", i, j, x);

  if(i < 0 || i >= a->dim[0] || j < 0 || j >= a->dim[1]) {
    PyErr_SetString(PyExc_IndexError, "indices out of range");
    return -1;
  }

  /* Find element to be set (or removed) */
  col = last = -1;
  k = a->root[i]; 
  while (k != -1) {
    col = a->col[k];
    if (col >= j)
      break;
    last = k;
    k = a->link[k];
  }

  if (x != 0.0) {
    
    if (col == j) {

      /* element already exist */
      a->val[k] = x;

    } else {
      /* new element */

      /* find location for new element */
      if (a->free != -1) {
	
	/* use element from the free chain */
	new_elem = a->free;
	a->free = a->link[new_elem];
	
      } else {
	
	/* append new element to the end */
	new_elem = a->nnz;
	
	/* test if there is space for a new element */
	if (a->nnz == a->nalloc) {
	  int nalloc_new;
	  
	  /* increase size of idx, val and link arrays */
	  nalloc_new = (int)((double)INCREASE_FACTOR * a->nalloc) + 1;
	  if ((temp = PyMem_Resize(a->col, int, nalloc_new)) == NULL)
	    goto fail;
	  else
	    a->col = temp;
	  if ((temp = PyMem_Resize(a->link, int, nalloc_new)) == NULL)
	    goto fail;
	  else
	    a->link = temp;
	  if ((temp = PyMem_Resize(a->val, double, nalloc_new)) == NULL)
	    goto fail;
	  else
	    a->val = temp;
	  a->nalloc = nalloc_new;
	}
      }

      a->val[new_elem] = x;
      a->col[new_elem] = j;
      a->link[new_elem] = k;
      if (last == -1)
	a->root[i] = new_elem;
      else
	a->link[last] = new_elem;
      
      a->nnz ++;
    }

  } else { /* x == 0.0 */
    
    if (col == j) {
      /* relink row i */
      if (last == -1)
	a->root[i] = a->link[k];
      else
	a->link[last] = a->link[k];
      /* add element to free list */
      a->link[k] = a->free;
      a->free = k;
      
      a->nnz --;
    }
  }
  return 0;
    
 fail:
  PyErr_NoMemory();
  return -1;
}

/*
 *  Update-add matrix entry: a[i,j] += x
 */

static int 
SpMatrix_LLMatUpdateItemAdd(LLMatObject *a, int i, int j, double x) {

  void *temp;
  int k, new_elem, col, last;

  if (a->issym && i < j) {
    PyErr_SetString(PyExc_IndexError, //SpMatrix_ErrorObject, 
		    "write operation to upper triangle of symmetric matrix");
    return -1;
  } 

  if (x == 0.0)
    return 0;

  /* Find element to be updated */
  col = last = -1;
  k = a->root[i]; 
  while (k != -1) {
    col = a->col[k];
    if (col >= j)
      break;
    last = k;
    k = a->link[k];    
  }
  if (col == j) {
    /* element already exists: compute updated value */
    x += a->val[k];

    if (x == 0.0) {
      /* the updated element is zero and must be removed */
    
      /* relink row i */
      if (last == -1)
	a->root[i] = a->link[k];
      else
	a->link[last] = a->link[k];
      /* add element to free list */
      a->link[k] = a->free;
      a->free = k;
	
      a->nnz --;
    } else {
      a->val[k] = x;
    }

  } else {
    /* new item */
    if (a->free != -1) {

      /* use element from the free chain */
      new_elem = a->free;
      a->free = a->link[new_elem];

    } else {

      /* append new element to the end */
      new_elem = a->nnz;

      /* test if there is space for a new element */
      if (a->nnz == a->nalloc) {
	int nalloc_new;

	/* increase size of idx, val and link arrays */
	nalloc_new = (int)((double)INCREASE_FACTOR * a->nalloc) + 1;
	if ((temp = PyMem_Resize(a->col, int, nalloc_new)) == NULL)
	  goto fail;
	else
	  a->col = temp;
	if ((temp = PyMem_Resize(a->link, int, nalloc_new)) == NULL)
	  goto fail;
	else
	  a->link = temp;
	if ((temp = PyMem_Resize(a->val, double, nalloc_new)) == NULL)
	  goto fail;
	else
	  a->val = temp;
	a->nalloc = nalloc_new;
      }
    }
  
    a->val[new_elem] = x;
    a->col[new_elem] = j;
    a->link[new_elem] = k;
    if (last == -1)
      a->root[i] = new_elem;
    else
      a->link[last] = new_elem;
    a->nnz ++;
  }
  return 0;
  
 fail:
  PyErr_NoMemory();
  return -1;
}

/*
 *  Shrink val, col and link arrays of matrix to minimal size
 */

static int
LLMat_Compress(LLMatObject *self, int *nofFreed) {

  double *val;
  int *col, *link;
  int i, k, k_next, k_last, k_new, nalloc_new;

  nalloc_new = self->nnz;   /* new size for val, col and link arrays */

  /* remove entries with k >= nalloc_new from free list */
  k_last = -1;
  k = self->free;
  while (k != -1) {
    k_next =  self->link[k];
    if (k >= nalloc_new) {
      if (k_last == -1)
	self->free = k_next;
      else
	self->link[k_last] = k_next;
    } else
      k_last = k;
    k = k_next;
  }
  /* reposition matrix entries with k >= nalloc_new */
  for (i = 0; i < self->dim[0]; i ++) {
    k_last = -1;
    for (k = self->root[i]; k != -1; k = self->link[k]) {
      if (k >= nalloc_new) {
	k_new = self->free;
	if (k_last == -1)
	  self->root[i] = k_new;
	else
	  self->link[k_last] = k_new;
	self->free = self->link[k_new];
	self->val[k_new] = self->val[k];
	self->col[k_new] = self->col[k];
	self->link[k_new] = self->link[k];
	k_last = k_new;
      } else 
	k_last = k;
    }
  }
  /* shrink arrays */
  if( !(col = PyMem_Resize(self->col, int, nalloc_new)) ) goto fail;
  self->col = col;
  if( !(link = PyMem_Resize(self->link, int, nalloc_new)) ) goto fail;
  self->link = link;
  if( !(val = PyMem_Resize(self->val, double, nalloc_new)) ) goto fail;
  self->val = val;

  *nofFreed = self->nalloc - nalloc_new;
  self->nalloc = nalloc_new;
  return 0;

 fail:
  PyErr_NoMemory();
  return -1;
}

/***********************************************************************/
/*  R o u t i n e s   f o r   h a n d l i n g   s u b m a t r i c e s  */
/***********************************************************************/

/* Create a index array from a Python list or slice */

long* create_indexlist(long *len, long maxlen, PyObject *A) {
  long *index;
  long  i, j;
  Py_ssize_t start, stop, step, length;
  PyObject *val;

  /* Integer */
  if( PyInt_Check(A) ) {
    i = PyInt_AS_LONG(A);
    if( (index = calloc(1, sizeof(long))) ) index[0] = i;
    *len = 1;
    return index;
  }

  /* Slice */
  if( PySlice_Check(A) ) {

    if( PySlice_GetIndicesEx((PySliceObject*)A, maxlen,
                             &start, &stop, &step, &length) < 0) return NULL;

    if( (index = calloc(length, sizeof(long))) )
      for( i=start, j=0; j<length; i += step, j++) index[j] = i;

    *len = (int)length;
    return index;
  }

  /* List */
  if( PyList_Check(A) ) {
    length = PyList_Size(A);
    if( !(index = calloc(length, sizeof(long))) ) return NULL;
    for( i=0; i<length; i++ ) {
      val = PyList_GetItem(A, (Py_ssize_t)i);
      if( PyInt_Check(val) )
        index[i] = PyInt_AS_LONG(val);
      else {
        free(index);
        PyErr_SetString(PyExc_ValueError, "Index must be a list of integers");
        return NULL;
      }
    }
    *len = (int)length;
    return index;
  }

  /* Numpy array */
  if( PyArray_Check(A) ) {
    npy_intp length0 = PyArray_DIM(A, 0);
    PyObject *iterator0 = PyArray_IterNew(A);

    if( !(index = calloc(length0, sizeof(long))) ) {
      Py_XDECREF(iterator0);
      return NULL;
    }

    PyArray_ITER_RESET(iterator0);
    i = 0;

    while( PyArray_ITER_NOTDONE(iterator0) ) {
      index[i] = *(long*)PyArray_ITER_DATA(iterator0);
      PyArray_ITER_NEXT(iterator0);
      i++;
    }
    *len = (int)length0;
    Py_DECREF(iterator0);
    return index;
  }

  PyErr_SetString(PyExc_TypeError, "Invalid index type");
  return NULL;
}

/*
 *  Copy submatrix src[start0:stop0,start1:stop1] into a newly allocated matrix.
 *  If transpose is non-zero, the transposed submatrix is copied.
 *
 *  Helper routine for get_submatrix
 */

static int copySubMatrix_FromList(LLMatObject *src, LLMatObject *dst,
                                  long *irow, int nrow, long *jcol, int ncol) {

  // Transfer elements (i,j) in irow X jcol from src to dst.
  // In Python notation: dst = src[irow,jcol].

  int i, j, row, col;
  double val;
  
  // Move elements over one by one
  for( i = 0; i < nrow; i++ ) {
    row = irow[i];
    for( j = 0; j < ncol; j++ ) {
      col = jcol[j];
      val = SpMatrix_LLMatGetItem(src, row, col);
      // Element (row,col) in src goes to location (i,j) in dst
      if( SpMatrix_LLMatSetItem(dst, i, j, val) ) {
        PyErr_SetString(PyExc_ValueError, "SpMatrix_LLMatSetItem failed");
        return -1;
      }
    }
  }
  return 0;
}

static PyObject *getSubMatrix_FromList(LLMatObject *self,
                                       PyObject *index0, PyObject *index1) {

  LLMatObject *dst;           // Destination matrix
  int dim[2];                 // Destination matrix shape
  int symmetric = 0;          // For now, always return a general matrix
  int res;

  // Both index sets are a single integer
  if( PyInt_Check(index0) && PyInt_Check(index1) ) {
    long row = PyInt_AS_LONG(index0);
    long col = PyInt_AS_LONG(index1);
    while( row < 0 ) row += self->dim[0];
    while( col < 0 ) col += self->dim[1];
    if( row >= self->dim[0] || col >= self->dim[1] ) {
      PyErr_SetString(PyExc_IndexError, "Index ouf of bounds");
      return NULL;
    }
    return PyFloat_FromDouble( SpMatrix_LLMatGetItem(self, row, col) );
  }

  // Both index sets are Python slices
  if( PySlice_Check(index0) && PySlice_Check(index1) ) {

    int i, j, k, row, col, t, inslice1;
    Py_ssize_t start0, stop0, step0, length0;
    Py_ssize_t start1, stop1, step1, length1;
    double val;

    if( PySlice_GetIndicesEx((PySliceObject*)index0, self->dim[0],
                             &start0, &stop0, &step0, &length0) < 0)
      return NULL;

    if( PySlice_GetIndicesEx((PySliceObject*)index1, self->dim[1],
                             &start1, &stop1, &step1, &length1) < 0)
      return NULL;

    dim[0] = length0; dim[1] = length1;
    dst = (LLMatObject *)SpMatrix_NewLLMatObject(dim, symmetric, self->nnz);
    if( !dst ) return NULL;

    row = start0;
    i = 0;    // Set row index of first element to be output

    // Scan each row in turn
    while( (step0 > 0 && row < stop0) || (step0 < 0 && row > stop0) ) {

      // Scan current row
      for( k = self->root[row]; k != -1; k = self->link[k] ) {

        col = self->col[k]; // Col index of current nonzero element
        inslice1 = 0;
        if( step1 > 0 && col >= start1 && col < stop1 )
          inslice1 = (col % step1 == 0);
        else if( step1 < 0 && col <= start1 && col > stop1 )
          inslice1 = (col % step1 == 0);

        if( inslice1 ) {    // Want this element

          if( self->issym && row < col ) {  // Swap row and col
            t = col; col = row; row = t;
          }

          val = SpMatrix_LLMatGetItem(self, row, col);

          j = (long)((col - start1)/step1); // Col index of element in output

          if( SpMatrix_LLMatSetItem(dst, i, j, val) ) {
            Py_DECREF(dst);
            PyErr_SetString(PyExc_ValueError, "SpMatrix_LLMatSetItem failed");
            return NULL;
          }

        }
      }

      row += step0;  // Move to next row
      i++;
    }

    /* Here is another version of slicing which is less efficient since
     * it scans the whole slices and extracts every single element in the
     * slices, including all the potentially zero elements. */
    /*
    row = start0;
    i = 0;
    while( (step0 > 0 && row < stop0) || (step0 < 0 && row > stop0) ) {
      col = start1;
      j = 0;
      while( (step1 > 0 && col < stop1) || (step1 < 0 && col > stop1) ) {
        val = SpMatrix_LLMatGetItem(self, row, col);
        if( SpMatrix_LLMatSetItem(dst, i, j, val) ) {
          Py_DECREF(dst);
          PyErr_SetString(PyExc_ValueError, "SpMatrix_LLMatSetItem failed");
          return NULL;
        }
        col += step1;
        j++;
      }
      row += step0;
      i++;
    }
    */

    /* A third version of slicing would be to build a column index for the
     * input matrix, scan both the row and column linked lists and retain
     * only common elements.
     */

    return (PyObject *)dst;
  }

  // Both index sets are Python lists
  if( PyList_Check(index0) && PyList_Check(index1) ) {

    int i, j;
    long row, col;
    Py_ssize_t length0 = PyList_Size(index0);
    Py_ssize_t length1 = PyList_Size(index1);
    PyObject *ind;
    double val;

    dim[0] = length0; dim[1] = length1;
    dst = (LLMatObject *)SpMatrix_NewLLMatObject(dim, symmetric, self->nnz);
    if( !dst ) return NULL;

    for( i = 0; i < length0; i++ ) {
      ind = PyList_GetItem(index0, (Py_ssize_t)i);
      if( PyInt_Check(ind) )
        row = PyInt_AS_LONG(ind);
      else {
        Py_DECREF(dst);
        PyErr_SetString(PyExc_ValueError, "Invalid list item");
        return NULL;
      }
      for( j = 0; j < length1; j++ ) {
        ind = PyList_GetItem(index1, (Py_ssize_t)j);
        if( PyInt_Check(ind) )
          col = PyInt_AS_LONG(ind);
        else {
          Py_DECREF(dst);
          PyErr_SetString(PyExc_ValueError, "Invalid list item");
          return NULL;
        }

        val = SpMatrix_LLMatGetItem(self, row, col);
        if( SpMatrix_LLMatSetItem(dst, i, j, val) ) {
          Py_DECREF(dst);
          PyErr_SetString(PyExc_ValueError, "SpMatrix_LLMatSetItem failed");
          return NULL;
        }
      }
    }

    return (PyObject *)dst;
  }

  // Both index sets are Numpy arrays
  if( PyArray_Check(index0) && PyArray_Check(index1) ) {

    int i, j;
    long row, col;
    double val;
    PyObject *iterator0 = PyArray_IterNew(index0);  // Use iterators because may
    PyObject *iterator1 = PyArray_IterNew(index1);  // not be contiguous
    npy_intp length0 = PyArray_DIM(index0, 0);
    npy_intp length1 = PyArray_DIM(index1, 0);

    dim[0] = length0; dim[1] = length1;
    dst = (LLMatObject *)SpMatrix_NewLLMatObject(dim, symmetric, self->nnz);
    if( !dst ) {
      Py_XDECREF(iterator0);
      Py_XDECREF(iterator1);
      return NULL;
    }

    PyArray_ITER_RESET(iterator0);
    i = 0;

    while( PyArray_ITER_NOTDONE(iterator0) ) {

      row = *(long*)PyArray_ITER_DATA(iterator0);

      PyArray_ITER_RESET(iterator1);
      j = 0;

      while( PyArray_ITER_NOTDONE(iterator1) ) {

        col = *(long*)PyArray_ITER_DATA(iterator1);
        val = SpMatrix_LLMatGetItem(self, row, col);

        if( SpMatrix_LLMatSetItem(dst, i, j, val) ) {
          Py_DECREF(dst);
          Py_DECREF(iterator0);
          Py_DECREF(iterator1);
          PyErr_SetString(PyExc_ValueError, "SpMatrix_LLMatSetItem failed");
          return NULL;
        }

        PyArray_ITER_NEXT(iterator1);
        j++;
      }

      PyArray_ITER_NEXT(iterator0);
      i++;
    }

    Py_DECREF(iterator0);
    Py_DECREF(iterator1);
    return (PyObject *)dst;
  }

  // If we have a mixture of index sets, gather both sets into arrays.
  {
      long *irow, *jcol;
      long  nrow,  ncol;

      // Create index list from first index
      if( !(irow = create_indexlist(&nrow, self->dim[0], index0)) ) {
        PyErr_SetString(PyExc_IndexError, "Error creating first index list");
        Py_INCREF(Py_None);
        return Py_None;
      }

      // Create index list from second index
      if( !(jcol = create_indexlist(&ncol, self->dim[1], index1)) ) {
        PyErr_SetString(PyExc_IndexError, "Error creating second index list");
        Py_INCREF(Py_None);
        return Py_None;
      }

      dim[0] = nrow; dim[1] = ncol;
      dst = (LLMatObject *)SpMatrix_NewLLMatObject(dim, symmetric, self->nnz);
      if( !dst ) return NULL;

      res = copySubMatrix_FromList(self, (LLMatObject*)dst, irow, nrow, jcol, ncol);
      if( res ) {
        Py_DECREF(dst);
        return NULL;
      }
      return (PyObject *)dst;
  }

}

/*
 *  Delete all non-zero entries from slice self[start0:stop0, start1:stop1]
 */

static int
clear_submatrix(LLMatObject *self,
                int start0, int stop0, int start1, int stop1) {

  int i, j, k, next, last;

  for (i = start0; i < stop0; i ++) {
    last = -1;
    k = self->root[i];
    while (k != -1) {
      j = self->col[k];
      next = self->link[k];
      if (start1 <= j && j < stop1) {
	/* remove element */
	if (last == -1)
	  self->root[i] = next;
	else
	  self->link[last] = next;
	/* add element to free list */
	self->link[k] = self->free;
	self->free = k;
	self->nnz--;
      } else
	last = k;
      k = next;
    }
  }
  return 0;
}

/*
 *  Assign slice self[irow,jcol] = other.
 */

//static int setSubMatrix_FromList(LLMatObject *self, LLMatObject *other,
//                                 long *irow, int nrow, long *jcol, int ncol) {
static void //int
setSubMatrix_FromList(LLMatObject *self, PyObject *other,
                                 PyObject *index0, PyObject *index1) {

  LLMatObject *mat = NULL;
  int other_is_num = 0, other_is_sym;
  double val = 0.0;

  if( PyInt_Check(other) ) {
    val = (double)PyInt_AsLong(other);
    other_is_num = 1;
  } else if( PyFloat_Check(other) ) {
    val = PyFloat_AsDouble(other);
    other_is_num = 1;
  }
  
  //other_is_num = PyArg_Parse(other, "d;array item must be float", &val);

  // Both index sets are a single integer
  if( PyInt_Check(index0) && PyInt_Check(index1) ) {

    long row = PyInt_AS_LONG(index0);
    long col = PyInt_AS_LONG(index1);

    if( !other_is_num ) { //return -1;
      PyErr_SetString(PyExc_ValueError, "Value must be double");
      return;
    }
    //return SpMatrix_LLMatSetItem(self, row, col, val);
    while( row < 0 ) row += self->dim[0];
    while( col < 0 ) col += self->dim[1];
    if( row >= self->dim[0] || col >= self->dim[1] ) {
      PyErr_SetString(PyExc_IndexError, "Index out of bounds");
      return;
    }
    SpMatrix_LLMatSetItem(self, row, col, val);
    return;
  }

  if( !other_is_num ) {
    mat = (LLMatObject *)other;
    other_is_sym = mat->issym;
  } else
    other_is_sym = 0;

  // Both index sets are Python slices
  if( PySlice_Check(index0) && PySlice_Check(index1) ) {

    long i, j, k, row, col, t;
    Py_ssize_t start0, stop0, step0, length0;
    Py_ssize_t start1, stop1, step1, length1;

    //printf("PySparse:: we have two slices...\n");

    if( PySlice_GetIndicesEx((PySliceObject*)index0, self->dim[0],
                             &start0, &stop0, &step0, &length0) < 0) {
      PyErr_SetString(PyExc_IndexError, "Erroneous indices");
      return;
      //return -1;
    }

    if( PySlice_GetIndicesEx((PySliceObject*)index1, self->dim[1],
                             &start1, &stop1, &step1, &length1) < 0) {
      PyErr_SetString(PyExc_IndexError, "Erroneous indices");
      return;
      //return -1;
    }

    //printf("New block of size (%ld,%ld)\n", length0, length1);

    if( other_is_num) {

      // Special case for A[slice,slice] = scalar
      row = start0;
      while( (step0 > 0 && row < stop0) || (step0 < 0 && row > stop0) ) {
        col = start1;
        while( (step1 > 0 && col < stop1) || (step1 < 0 && col > stop1) ) {
          if( self->issym && row < col ) {
            PyErr_SetString(PyExc_IndexError,  //SpMatrix_ErrorObject,
                            "Writing to upper triangle of symmetric matrix");
            return; // -1;
          }
          if( SpMatrix_LLMatSetItem(self, row, col, val) ) {
            PyErr_SetString(PyExc_ValueError, "SpMatrix_LLMatSetItem failed");
            return; // -1;
          }
          col += step1;
        }
        row += step0;
      }

    } else {

      if( mat->dim[0] != length0 || mat->dim[1] != length1 ) {
        printf("In LHS=RHS, RHS has shape (%d,%d), LHS has shape (%d,%d)\n",
               mat->dim[0], mat->dim[1], (int)length0, (int)length1);
        PyErr_SetString(PyExc_ValueError, "Matrix shapes are different");
        return; // -1;
      }

      row = start0;
      i = 0;    // Set row index of first element to be assigned

      // Scan each row in turn
      while( (step0 > 0 && row < stop0) || (step0 < 0 && row > stop0) ) {

        // Scan current row in matrix to be assigned
        for( k = mat->root[i]; k != -1; k = mat->link[k] ) {
        
          j = mat->col[k]; // Col index of current nonzero element
          col = start1 + j * step1; // Col index to be assigned to

          if( mat->issym && row < col ) {  // Swap row and col
            t = col; col = row; row = t;
          }

          // Ensure write operation is permitted
          if( self->issym && row < col ) {
            PyErr_SetString(PyExc_IndexError, //SpMatrix_ErrorObject,
                            "Writing to upper triangle of symmetric matrix");
            return; // -1;
          }
          
          val = SpMatrix_LLMatGetItem(mat, i, j);
          //printf("  (%ld,%ld) -> (%ld,%ld). Val = %g\n", i,j,row,col,val);

          if( SpMatrix_LLMatSetItem(self, row, col, val) ) {
            PyErr_SetString(PyExc_ValueError, "SpMatrix_LLMatSetItem failed");
            return; // -1;
          }
        }

        row += step0;  // Move to next row
        i++;
      }
    }
    return; // 0;
  }

  // Both index sets are Python lists
  if( PyList_Check(index0) && PyList_Check(index1) ) {

    //printf("PySparse:: We have two lists...\n");

    long i, j, row, col, t;
    Py_ssize_t length0 = PyList_Size(index0);
    Py_ssize_t length1 = PyList_Size(index1);
    PyObject *ind;

    if( !other_is_num )
      if( mat->dim[0] != length0 || mat->dim[1] != length1 ) {
        PyErr_SetString(PyExc_ValueError, "Matrix shapes are different");
        return; // -1;
      }

    for( i = 0; i < length0; i++ ) {
      ind = PyList_GetItem(index0, (Py_ssize_t)i);
      if( PyInt_Check(ind) )
        row = PyInt_AS_LONG(ind);
      else {
        PyErr_SetString(PyExc_IndexError, "Invalid list item");
        return; // -1;
      }
      for( j = 0; j < length1; j++ ) {
        ind = PyList_GetItem(index1, (Py_ssize_t)j);
        if( PyInt_Check(ind) )
          col = PyInt_AS_LONG(ind);
        else {
          PyErr_SetString(PyExc_IndexError, "Invalid list item");
          return; // -1;
        }

        if( !other_is_num )
          val = SpMatrix_LLMatGetItem(mat, i, j);

        if( other_is_sym && row < col ) {
          t = col;
          col = row;
          row = t;
        }

        //printf("Val = %g goes to position (%ld,%ld)\n", val, row, col);

        // Ensure write operation is permitted
        if( self->issym && row < col ) {
          PyErr_SetString(PyExc_IndexError, //SpMatrix_ErrorObject,
                          "Writing to upper triangle of symmetric matrix");
          return; // -1;
        }

        if( SpMatrix_LLMatSetItem(self, row, col, val) ) {
          PyErr_SetString(PyExc_ValueError, "SpMatrix_LLMatSetItem failed");
          return; // -1;
        }
      }
    }

    return; // 0;
  }

  // Both index sets are Numpy arrays
  if( PyArray_Check(index0) && PyArray_Check(index1) ) {

    long i, j, row, col, t;
    PyObject *iterator0 = PyArray_IterNew(index0);  // Use iterators because may
    PyObject *iterator1 = PyArray_IterNew(index1);  // not be contiguous
    npy_intp length0 = PyArray_DIM(index0, 0);
    npy_intp length1 = PyArray_DIM(index1, 0);

    if( !other_is_num )
      if( mat->dim[0] != length0 || mat->dim[1] != length1 ) {
        Py_DECREF(iterator0);
        Py_DECREF(iterator1);
        PyErr_SetString(PyExc_ValueError, "Matrix shapes are different");
        return; // -1;
      }

    PyArray_ITER_RESET(iterator0);
    i = 0;

    while( PyArray_ITER_NOTDONE(iterator0) ) {

      row = *(long*)PyArray_ITER_DATA(iterator0);

      PyArray_ITER_RESET(iterator1);
      j = 0;

      while( PyArray_ITER_NOTDONE(iterator1) ) {

        col = *(long*)PyArray_ITER_DATA(iterator1);

        if( !other_is_num )
          val = SpMatrix_LLMatGetItem(mat, i, j);

        if( other_is_sym && row < col ) {
          t = col;
          col = row;
          row = t;
        }

        // Ensure write operation is permitted
        if( self->issym && row < col ) {
          Py_DECREF(iterator0);
          Py_DECREF(iterator1);
          PyErr_SetString(PyExc_IndexError, //SpMatrix_ErrorObject,
                          "Writing to upper triangle of symmetric matrix");
          return; // -1;
        }

        if( SpMatrix_LLMatSetItem(self, row, col, val) ) {
          Py_DECREF(iterator0);
          Py_DECREF(iterator1);
          PyErr_SetString(PyExc_ValueError, "SpMatrix_LLMatSetItem failed");
          return; // -1;
        }

        PyArray_ITER_NEXT(iterator1);
        j++;
      }

      PyArray_ITER_NEXT(iterator0);
      i++;
    }

    Py_DECREF(iterator0);
    Py_DECREF(iterator1);
    return; // 0;
  }

  // If we have a mixture of index sets, gather both sets into arrays.
  {
    long *irow, *jcol;
    long  i, j, nrow,  ncol, row, col;

    // Create index list from first index
    if( !(irow = create_indexlist(&nrow, self->dim[0], index0)) ) {
      PyErr_SetString(PyExc_IndexError, "Error creating first index list");
      return; // -1;
    }

    // Create index list from second index
    if( !(jcol = create_indexlist(&ncol, self->dim[1], index1)) ) {
      PyErr_SetString(PyExc_IndexError, "Error creating second index list");
      return; // -1;
    }

    if( !other_is_num )
      if( mat->dim[0] != nrow || mat->dim[1] != ncol ) {
        PyErr_SetString(PyExc_ValueError, "Matrix shapes are different");
        return; // -1;
      }

    for( i = 0; i < nrow; i++ ) {
      row = irow[i];
      for( j = 0; j < ncol; j++ ) {
        col = jcol[j];

        if( other_is_sym && row < col ) { 
          col = row;
          row = jcol[j];
        }

        // Ensure write operation is permitted
        if( self->issym && row < col ) {
          PyErr_SetString(PyExc_IndexError, //SpMatrix_ErrorObject,
                          "Writing to upper triangle of symmetric matrix");
          return; // -1;
        }

        // Insert element into self
        if( !other_is_num )
          val = SpMatrix_LLMatGetItem(mat, i, j);

        if( SpMatrix_LLMatSetItem(self, row, col, val) ) {
          PyErr_SetString(PyExc_ValueError, "SpMatrix_LLMatSetItem failed");
          return; // -1;
        }
      }
    }
    return; // 0;
  }
}

/*****************************************************************************/
/*  M a t r i x - v e c t o r   m u l t i p l i c a t i o n   k e r n e l s  */
/*****************************************************************************/

static void
ll_matvec_kernel(int m, double *x, double *y,
		 double *val, int *col, int *link, int *root) {
  double s;
  int i, k;
  
  for (i = 0; i < m; i ++) {
    s = 0.0;
    k = root[i];
    while (k != -1) {
      s += val[k] * x[col[k]];
      k = link[k];
    }
    y[i] = s;
  }
}

static void
ll_matvec_kernel_stride(int m, 
			double *x, int incx, 
			double *y, int incy,
			double *val, int *col, int *link, int *root) {
  double s;
  int i, k;
  
  for (i = 0; i < m; i ++) {
    s = 0.0;
    k = root[i];
    while (k != -1) {
      s += val[k] * x[col[k]*incx];
      k = link[k];
    }
    y[i*incy] = s;
  }
}

static void
ll_matvec_kernel_sym(int m, double *x, double *y,
		     double *val, int *col, int *link, int *root) {
  double s, v, xi;
  int i, j, k;
  
  for (i = 0; i < m; i ++) {
    xi = x[i];
    s = 0.0;
    k = root[i];
    while (k != -1) {
      j = col[k];
      v = val[k];
      s += v * x[j];
      if (i != j)
	y[j] += v * xi;
      k = link[k];
    }
    y[i] = s;
  }
}

static void
ll_matvec_kernel_stride_sym(int m, 
			    double *x, int incx, 
			    double *y, int incy,
			    double *val, int *col, int *link, int *root) {
  double s, v, xi;
  int i, j, k;
  
  for (i = 0; i < m; i ++) {
    xi = x[i*incx];
    s = 0.0;
    k = root[i];
    while (k != -1) {
      j = col[k];
      v = val[k];
      s += v * x[j*incx];
      if (i != j)
	y[j*incy] += v * xi;
      k = link[k];
    }
    y[i*incy] = s;
  }
}

static void
ll_matvec_transp_kernel(int m, int n, double *x, double *y,
			    double *val, int *col, int *link, int *root) {
  double xi;
  int i, k;
  
  for (i = 0; i < n; i ++)
    y[i] = 0.0;
  
  for (i = 0; i < m; i ++) {
    xi = x[i];
    k = root[i];
    while (k != -1) {
      y[col[k]] += val[k] * xi;
      k = link[k];
    }
  }
}

static void
ll_matvec_transp_kernel_stride(int m, int n, 
			       double *x, int incx, 
			       double *y, int incy,
			       double *val, int *col, int *link, int *root) {
  double xi;
  int i, k;
  
  for (i = 0; i < n; i ++)
    y[i*incy] = 0.0;
  
  for (i = 0; i < m; i ++) {
    xi = x[i*incx];
    k = root[i];
    while (k != -1) {
      y[col[k]*incy] += val[k] * xi;
      k = link[k];
    }
  }
}

/*********************************************/
/*  L L M a t   o b j e c t   m e t h o d s  */
/*********************************************/

static char LLMat_matvec_transp_doc[] = "A.matvec_transp(x, y)\n\
\n\
compute the sparse matrix-vector product y := A^T * x. \n\
A^T is the transpose of A, which is a d1 by d2 sparse matrix.\n\
x and y are two 1-dimensional Numpy arrays of appropriate size.";

static PyObject *
LLMat_matvec_transp(LLMatObject *self, PyObject *args)
{
  PyArrayObject *xp, *yp;
  size_t sd = sizeof(double);

  SPMATRIX_PARSE_ARGS_ARR_ARR_STRIDE(args, xp, yp, self->dim[0], self->dim[1]);
  
  if (xp->flags & CONTIGUOUS &&  yp->flags & CONTIGUOUS)
    if (self->issym)
      ll_matvec_kernel_sym(self->dim[0],
                           (double *)(xp->data), (double *)(yp->data), 
			   self->val, self->col, self->link, self->root);
    else
      ll_matvec_transp_kernel(self->dim[0], self->dim[1],
                              (double *)(xp->data), (double *)(yp->data), 
			      self->val, self->col, self->link, self->root);
  else
    if (self->issym)
      ll_matvec_kernel_stride_sym(self->dim[0], 
				  (double *)(xp->data), xp->strides[0] / sd,
				  (double *)(yp->data), yp->strides[0] / sd,
				  self->val, self->col, self->link, self->root);
    else
      ll_matvec_transp_kernel_stride(self->dim[0], self->dim[1],
				     (double *)(xp->data), xp->strides[0] / sd,
				     (double *)(yp->data), yp->strides[0] / sd,
				     self->val, self->col, self->link,
                                     self->root);

  Py_INCREF(Py_None); 
  return Py_None;
}

static char LLMat_matvec_doc[] = "A.matvec(x, y)\n\
\n\
compute the sparse matrix-vector product y := A * x. \n\
A is a d1 by d2 sparse matrix.\n\
x and y are two 1-dimensional Numpy arrays of appropriate size.";

static PyObject *
LLMat_matvec(LLMatObject *self, PyObject *args)
{
  PyArrayObject *xp, *yp;
  size_t sd = sizeof(double);

  SPMATRIX_PARSE_ARGS_ARR_ARR_STRIDE(args, xp, yp, self->dim[1], self->dim[0]);
     
  if (xp->flags & CONTIGUOUS &&  yp->flags & CONTIGUOUS)
    if (self->issym)
      ll_matvec_kernel_sym(self->dim[0],
                           (double *)(xp->data), (double *)(yp->data), 
			   self->val, self->col, self->link, self->root);
    else
      ll_matvec_kernel(self->dim[0],
                       (double *)(xp->data), (double *)(yp->data), 
		       self->val, self->col, self->link, self->root);
  else
    if (self->issym)
      ll_matvec_kernel_stride_sym(self->dim[0], 
				  (double *)(xp->data), xp->strides[0] / sd,
				  (double *)(yp->data), yp->strides[0] / sd,
				  self->val, self->col, self->link, self->root);
    else
      ll_matvec_kernel_stride(self->dim[0], 
			      (double *)(xp->data), xp->strides[0] / sd,
			      (double *)(yp->data), yp->strides[0] / sd,
			      self->val, self->col, self->link, self->root);

  Py_INCREF(Py_None); 
  return Py_None;
}

/* Apply in-place column scaling to a matrix of type LL_Mat */

static char col_scale_doc[] = "A.col_scale(v)\n\
\n\
Scale the i-th column of A by v[i] in place for i=0, ..., ncol-1.";

static PyObject *
LLMat_col_scale(LLMatObject *self, PyObject *args) {

  PyObject *vIn;
  PyArrayObject *v = NULL;
  struct llColIndex *colIdx;
  double val;
  int j, k;

  // Read scale vector v.
  if( !PyArg_ParseTuple(args, "O", &vIn) ) {
    PyErr_SetString(SpMatrix_ErrorObject, "Cannot read input vector.");
    return NULL;
  }
  v = (PyArrayObject *)PyArray_ContiguousFromObject(vIn, PyArray_DOUBLE, 1, 1);
  if( v == NULL ) {
    PyErr_SetString(SpMatrix_ErrorObject, "Supply scaling vector as input.");
    return NULL;
  }

  // Check for dimensions mismatch.
  if( v->dimensions[0] != self->dim[1] ) {
    PyErr_SetString(SpMatrix_ErrorObject,
                    "Column scaling vector has wrong dimension.");
    return NULL;
  }

  // Build column index.
  if( SpMatrix_LLMatBuildColIndex(&colIdx, self, 1) ) {
    PyErr_SetString(SpMatrix_ErrorObject, "Cannot build column index.");
    return NULL;
  }

  // Process each column in turn
  for( j = 0; j < self->dim[1] ; j++ ) {
    val = ((double *)v->data)[j];

    // Scan column j.
    k = colIdx->root[j];
    while( k != -1 ) {
      self->val[k] *= val;
      k = colIdx->link[k];
    }
  }

  SpMatrix_LLMatDestroyColIndex(&colIdx);

  Py_INCREF(Py_None);
  return Py_None;
}

/* Apply in-place row scaling to a matrix of type LL_Mat */

static char row_scale_doc[] = "A.row_scale(v)\n\
\n\
Scale the i-th row of A by v[i] in place for i=0, ..., nrow-1.";

static PyObject *
LLMat_row_scale(LLMatObject *self, PyObject *args) {

  double val;
  int i, k;
  PyObject *vIn;
  PyArrayObject *v = NULL;

  // Read scale vector v.
  if( !PyArg_ParseTuple(args, "O", &vIn) ) {
    PyErr_SetString(SpMatrix_ErrorObject, "Cannot read input vector.");
    return NULL;
  }
  v = (PyArrayObject *)PyArray_ContiguousFromObject(vIn, PyArray_DOUBLE, 1, 1);
  if( v == NULL ) {
    PyErr_SetString(SpMatrix_ErrorObject, "Supply scaling vector as input.");
    return NULL;
  }

  // Check for dimensions mismatch.
  if( v->dimensions[0] != self->dim[0] ) {
    PyErr_SetString(SpMatrix_ErrorObject,
                    "Row scaling vector has wrong dimension.");
    return NULL;
  }

  // Process each row in turn.
  for( i = 0; i < self->dim[0]; i++ ) {
    val = ((double *)v->data)[i];

    // Scan row i.
    k = self->root[i];
    while( k != -1 ) {
      self->val[k] *= val;
      k = self->link[k];
    }
  }

  Py_INCREF(Py_None);
  return Py_None;
}


static char to_csr_doc[] = "A.to_csr()\n\
\n\
return a new CSRMatObject constructed from data of A";

static PyObject *
LLMat_to_csr(LLMatObject *self, PyObject *args)
{
  CSRMatObject *op;
  int i, j, k, r;
  
  if (!PyArg_ParseTuple(args, ""))
    return NULL;

  if (self->issym) {  /* Symmetric case */
    struct llColIndex *colIdx;

    if (SpMatrix_LLMatBuildColIndex(&colIdx, self, 0)) return NULL;
    assert(colIdx->nzUp == 0);

    op = (CSRMatObject *)newCSRMatObject(self->dim,
                                         2*colIdx->nzLo + colIdx->nzDiag);
    if (op == NULL) {
      SpMatrix_LLMatDestroyColIndex(&colIdx);
      return NULL;
    }

    r = 0;
    op->ind[0] = 0;
    for (i = 0; i < self->dim[0]; i ++) {

      /* store self[0:i+1,i] in op[0:i+1,i] */
      k = self->root[i];
      while (k != -1) {
	op->val[r] = self->val[k];
	op->col[r] = self->col[k];
	r ++;
	k = self->link[k];
      }

      /* store self[i,i+1:n] in op[i+1:n,i] */
      k = colIdx->root[i];
      while (k != -1) {
	j = colIdx->row[k];
	op->val[r] = self->val[k];
	op->col[r] = j;
	r ++;
	k = colIdx->link[k];
      }
      
      op->ind[i+1] = r;
    }
  
    SpMatrix_LLMatDestroyColIndex(&colIdx);

  } else {  /* Unsymmetric case */

    op = (CSRMatObject *)newCSRMatObject(self->dim, self->nnz);
    if (op == NULL)
      return NULL;

    r = 0;
    op->ind[0] = 0;
    for (i = 0; i < self->dim[0]; i ++) {
      k = self->root[i];
      while (k != -1) {
	op->val[r] = self->val[k];
	op->col[r] = self->col[k];
	r ++;
	k = self->link[k];
      }
      op->ind[i+1] = r;
    }
  }

  return (PyObject *)op;
}

static char to_sss_doc[] = "a.to_sss()\n\
\n\
return a new SSSMatObject constructed from the lower triangle of a";

static PyObject *
LLMat_to_sss(LLMatObject *self, PyObject *args)
{
  SSSMatObject *op;
  int i, j, k, r, n, nnz;
  
  if (!PyArg_ParseTuple(args, ""))
    return NULL;

  /* test for square matrix */
  n = self->dim[0];
  if (n != self->dim[1]) {
    PyErr_SetString(PyExc_ValueError, "Matrix must be square");
    return NULL;
  }
  
  /* 1st pass: compute number of non-zeros in lower triangle */
  nnz = 0;
  for (i = 0; i < n; i ++) {
    k = self->root[i];
    while (k != -1) {
      j = self->col[k];
      if (i > j)
	nnz ++;
      k = self->link[k];
    }
  }
  
  /* allocate new SSS matrix */
  op = (SSSMatObject *)newSSSMatObject(n, nnz);
  if (op == NULL)
    return NULL;

  /* 2nd pass: fill SSSMatObject */
  for (i = 0; i < n; i ++)
    op->diag[i] = 0.0;
  r = 0;
  op->ind[0] = 0;
  for (i = 0; i < n; i ++) {
    k = self->root[i];
    while (k != -1) {
      j = self->col[k];
      if (i > j) {
	op->val[r] = self->val[k];
	op->col[r] = j;
	r ++;
      } else if (i == j)
	op->diag[i] = self->val[k];
      k = self->link[k];
    }
    op->ind[i+1] = r;
  }

  return (PyObject *)op;
}

static char LLMat_generalize_doc[] = \
"convert ll_mat object from symmetric to non-symmetric form (in-place).";

static PyObject *
LLMat_generalize(LLMatObject *self, PyObject *args) {
  int i, j, k;

  if (!PyArg_ParseTuple(args, "")) return NULL;

  if (self->issym) {
    self->issym = 0;
    for (i = 0; i < self->dim[0]; i ++) {
      /* New elements are inserted into the matrix while it is being traversed.
	 However, this should not be a problem */
      for (k = self->root[i]; k != -1; k = self->link[k]) {
	j = self->col[k];
	if (i > j)
	  if (SpMatrix_LLMatSetItem(self, j, i, self->val[k]))
	    return NULL;
      }
    }
  }

  Py_INCREF(Py_None); 
  return Py_None;
}

static char LLMat_compress_doc[] = "A.compress() frees memory by reclaiming unused space in the data structures of A. Returns number of elements freed.";

static PyObject *
LLMat_compress(LLMatObject *self, PyObject *args) {
  int nofFreed;

  if (!PyArg_ParseTuple(args, "")) return NULL;
  if (LLMat_Compress(self, &nofFreed)) return NULL;
  return PyInt_FromLong(nofFreed);
}

static char export_mtx_doc[] = "A.export_mtx(fileName, precision=16)\n\
\n\
write matrix in Matrix Market format to fileName.\n\
\n\
Parameters:\n\
\n\
fileName:  string, name of the file to be created.\n\
precision: number of significant digits to be written for double values.";

static PyObject *
LLMat_export_mtx(LLMatObject *self, PyObject *args) {

  char *fileName;
  int precision = 16;
  MM_typecode matcode;
  FILE *f;
  int ret, i, k;

  if (!PyArg_ParseTuple(args, "s|i", &fileName, &precision)) return NULL;

  if( !(f = fopen(fileName, "w")) ) return PyErr_SetFromErrno(PyExc_IOError);
  
  mm_set_matrix(matcode); mm_set_sparse(matcode); mm_set_real(matcode);
  self->issym ? mm_set_symmetric(matcode) : mm_set_general(matcode);
    
  ret = mm_write_banner(f, matcode);
  if (ret) {
    PyErr_SetString(SpMatrix_ErrorObject, "Error writing file header");    
    return NULL;
  }
  
  ret = fprintf(f, "%% file created by pysparse module\n");
  if (ret < 0) {
    PyErr_SetString(PyExc_IOError, "Error writing file header");    
    return NULL;
  }

  ret = mm_write_mtx_crd_size(f, self->dim[0], self->dim[1], self->nnz);
  if (ret) {
    PyErr_SetString(SpMatrix_ErrorObject, "Error writing file header");    
    return NULL;
  }
  
  for (i = 0; i < self->dim[0]; i ++) {
    k = self->root[i];
    while (k != -1) {
      ret = fprintf(f, "%d %d %.*e\n",
                    i+1, self->col[k]+1, precision-1, self->val[k]);
      if (ret < 0) {
	PyErr_SetString(PyExc_IOError, "Error writing matrix data");    
	return NULL;
      }
      k = self->link[k];
    }
  }
  
  ret = fclose(f);
  if (ret)
    return PyErr_SetFromErrno(PyExc_IOError);

  Py_INCREF(Py_None); 
  return Py_None;  
}

static char copy_doc[] = "A.copy()\n\
\n\
return a (deep) copy of the matrix A.";

static PyObject *
LLMat_copy(LLMatObject *self, PyObject *args) {

  LLMatObject *new;
  int i, k;

  if (!PyArg_ParseTuple(args, "")) return NULL;
  
  new = (LLMatObject *)SpMatrix_NewLLMatObject(self->dim,
                                               self->issym, self->nnz);
  if (new == NULL) return NULL;

  for (i = 0; i < self->dim[0]; i++) {
    k = self->root[i];
    while (k != -1) {
      if (SpMatrix_LLMatSetItem(new, i, self->col[k], self->val[k]) == -1) {
	Py_DECREF(new);
	return NULL;
      }
      k = self->link[k];
    }
  }
  return (PyObject *)new;
}

static char update_add_at_doc[] = "A.update_add_at(val,id1,id2)\n\
\n\
for i in range(len(val)):\n\
    A[id1[i],id2[i]] += val[i]";

static PyObject *
LLMat_update_add_at(LLMatObject *self, PyObject *args) {
  PyObject *bIn;
  PyObject *id1in;
  PyObject *id2in;
  PyArrayObject *b = NULL;
  PyArrayObject *id1 = NULL;
  PyArrayObject *id2 = NULL;
  double v;
  int lenb, i, err;

  if (!PyArg_ParseTuple(args, "OOO", &bIn, &id1in, &id2in)) return NULL;

  b = (PyArrayObject *)PyArray_ContiguousFromObject(bIn, PyArray_DOUBLE, 1,1);
  if (b == NULL) goto fail;

  lenb = b->dimensions[0];

  id1 = (PyArrayObject *)PyArray_ContiguousFromObject(id1in, PyArray_LONG, 1,1);
  if (id1 == NULL) goto fail;
  
  id2 = (PyArrayObject *)PyArray_ContiguousFromObject(id2in, PyArray_LONG, 1,1);
  if (id2 == NULL) goto fail;

  if (lenb < 0 ) {
    PyErr_SetString(PyExc_IndexError, "Vector b has a negative size");
    goto fail;
  }

  if (id1->dimensions[0] != lenb ) {
    PyErr_SetString(PyExc_IndexError, "id1 does not have the same size as b");
    goto fail;
  }

  if (id2->dimensions[0] != lenb ) {
    PyErr_SetString(PyExc_IndexError, "id2 does not have the same size as b");
    goto fail;
  }

  /* Perform update add operation */
  for (i = 0; i < lenb; i ++) {
    v = ((double *)b->data)[i];
    err = SpMatrix_LLMatUpdateItemAdd(self, ((long *)id1->data)[i],
                                      ((long *)id2->data)[i], v);
    if( err == -1 ) goto fail;
  }

  Py_XDECREF(b);
  Py_XDECREF(id1);
  Py_XDECREF(id2);
  Py_INCREF(Py_None); 
  return Py_None;

 fail:
  if(b)   { Py_XDECREF(b);   }
  if(id1) { Py_XDECREF(id1); }
  if(id2) { Py_XDECREF(id2); }
  return NULL;
}

static char LLMat_norm_doc[] = "return p-norm of matrix\n\
\n\
A.norm(p) returns the p-norm of matrix A\n\
\n\
p is a string identifying the type of norm to be returned\n\
\n\
  '1'   -- return the 1-norm of A\n\
  'inf' -- return the infinity norm of A\n\
  'fro' -- return the frobenius norm of A";

static PyObject *
LLMat_norm(LLMatObject *self, PyObject *args)
{
  char *normType;

  struct llColIndex *colIdx;
  double norm, s, v;
  int i, k;
  
  if (!PyArg_ParseTuple(args, "s", &normType))
    return NULL;
  
  if (strcmp(normType, "1") == 0) {   /* l1 norm */

    if (self->issym) {
      PyErr_SetString(PyExc_NotImplementedError,
                      "Not implemented for symmetric matrices");
      return NULL;
    } else {

      if (SpMatrix_LLMatBuildColIndex(&colIdx, self, 1)) 
        return NULL;
      for( norm = 0.0, i = 0; i < self->dim[1]; i++ ) {
	for( s = 0.0, k = colIdx->root[i]; k != -1; k = colIdx->link[k] )
	  s += fabs(self->val[k]);
        norm = s > norm ? s : norm;
      }
      SpMatrix_LLMatDestroyColIndex(&colIdx);
    }

  } else if (strcmp(normType, "inf") == 0) {  /* infinity norm */

    if (self->issym) {
      PyErr_SetString(PyExc_NotImplementedError,
                      "Not implemented for symmetric matrices");
      return NULL;
    } else {
      for( norm = 0.0, i = 0; i < self->dim[0]; i++ ) {
	for( s = 0.0, k = self->root[i]; k != -1; k = self->link[k] )
	  s += fabs(self->val[k]);
        norm = s > norm ? s : norm;
      }
    }

  } else if (strcmp(normType, "fro") == 0) {  /* Frobenius norm */

    for( norm = 0.0, i = 0; i < self->dim[0]; i++ )
      for (k = self->root[i]; k != -1; k = self->link[k]) {
	v = self->val[k];
	norm += v*v;
	if (self->issym && self->col[k] != i) norm += v*v;
      }
    norm = sqrt(norm);

  } else {
    PyErr_SetString(PyExc_ValueError, "unknown norm type");
    return NULL;
  }

  return Py_BuildValue("d", norm);
}

static char shift_doc[] = "A.shift(sigma, B)\n\
\n\
shift the matrix:\n\
compute A = A + sigma * B\n\
where sigma is a scalar and B is a matrix of compatible size.";

static PyObject *
LLMat_shift(LLMatObject *self, PyObject *args) {
  LLMatObject *mat;
  double sigma, v;
  int i, j, k, err;
  
  if (!PyArg_ParseTuple(args, "dO!", &sigma, &LLMatType, &mat))
    return NULL;
  if (self->dim[0] != mat->dim[0] || self->dim[1] != mat->dim[1]) {
    PyErr_SetString(PyExc_ValueError, "matrix shapes do not match");
    return NULL;
  }

  if (self->issym == mat->issym) {
    for (i = 0; i < mat->dim[0]; i ++) {
      k = mat->root[i];
      while (k != -1) {
	err = SpMatrix_LLMatUpdateItemAdd(self, i,
                                          mat->col[k], sigma * mat->val[k]);
        if( err == -1 ) return NULL;
	k = mat->link[k];
      }
    }
  } else if (mat->issym) {
    for (i = 0; i < mat->dim[0]; i ++) {
      k = mat->root[i];
      while (k != -1) {
	j = mat->col[k];
	v = sigma * mat->val[k];
	if (SpMatrix_LLMatUpdateItemAdd(self, i, j, v) == -1)
	  return NULL;
	if (i != j)
	  if (SpMatrix_LLMatUpdateItemAdd(self, j, i, v) == -1)
	    return NULL;
	
	k = mat->link[k];
      }
    }
  } else {
    PyErr_SetString(PyExc_NotImplementedError, 
		    "Cannot shift symmetric matrix by non-symmetric matrix.");
    return NULL;
  }

  Py_INCREF(Py_None); 
  return Py_None;
}


static char keys_doc[] = "A.keys()\n\
\n\
Return a list of tuples (i,j) of non-zero matrix entries.";

/* static PyObject * */
/* LLMat_keys(LLMatObject *a, PyObject *args) { */
/*   PyObject *list;             /\* the list that will hold the keys *\/ */
/*   int i, j, k; */
/*   int pos = 0;                /\* position in list *\/ */
    
/*   if (!PyArg_ParseTuple(args, "")) return NULL; */
    
/*   if (!a->issym) { */

/*     if ((list = PyList_New(a->nnz)) == NULL) return NULL; */
        
/*     for (i = 0; i < a->dim[0]; i ++) { */
/*       k = a->root[i]; */
/*       while (k != -1) { */
/*         j = a->col[k]; */
/*                 PyList_SET_ITEM(list, pos++, Py_BuildValue("ii", i, j)); */
/*                 k = a->link[k]; */
/* 	    } */
/*         } */
/*         return list; */
        
/*     } else { */
/*         PyErr_SetString(PyExc_NotImplementedError,  */
/*                         "keys() doesn't yet support symmetric matrices"); */
/*         return NULL; */
/*     } */
/* } */

static PyObject *
LLMat_keys(LLMatObject *a, PyObject *args) {
  PyObject *listi;             /* the list that will hold keys i */
  PyObject *listj;             /* the list that will hold keys j */
  PyObject *list;             /* the list that will hold the keys */
  int i, j, k;
  int pos = 0;                /* position in list */
    
  if (!PyArg_ParseTuple(args, "")) return NULL;
    
  if (!a->issym) {

    //    printf("nnz is %i\n", a->nnz);
    if ((list = PyList_New(2)) == NULL) return NULL;
    if ((listi = PyList_New(a->nnz)) == NULL) return NULL;
    if ((listj = PyList_New(a->nnz)) == NULL) return NULL;
        
    for (i = 0; i < a->dim[0]; i ++) {
      k = a->root[i];
      while (k != -1) {
        j = a->col[k];
	PyList_SET_ITEM(listi, pos, PyInt_FromLong(i));
	PyList_SET_ITEM(listj, pos, PyInt_FromLong(j));
	pos++;
	k = a->link[k];
      }
    }

    PyList_SET_ITEM(list, 0, listi);
    PyList_SET_ITEM(list, 1, listj);

    return list;
        
  } else {
    PyErr_SetString(PyExc_NotImplementedError, 
		    "keys() doesn't yet support symmetric matrices");
    return NULL;
  }
}

static char values_doc[] = "A.values()\n\
\n\
Return a list of the non-zero matrix entries as floats.";

static PyObject *
LLMat_values(LLMatObject *a, PyObject *args) {
    PyObject *list;           /* the list that will hold the values */
    int i, k;
    int pos = 0;                /* position in list */

    if (!PyArg_ParseTuple(args, ""))
        return NULL;
    
    if (!a->issym) {
        
        if ((list = PyList_New(a->nnz)) == NULL)
            return NULL;
        
        for (i = 0; i < a->dim[0]; i ++) {
	    k = a->root[i];
	    while (k != -1) {
                PyList_SET_ITEM(list, pos++, PyFloat_FromDouble(a->val[k]));
                k = a->link[k];
	    }
        }
        return list;

    } else {
        PyErr_SetString(PyExc_NotImplementedError, 
                        "values() doesn't yet support symmetric matrices");
        return NULL;
    }
}

static char items_doc[] = "A.items()\n\
\n\
Return a list of tuples (indices, value) of\n\
the non-zero matrix entries' keys and values.\n\
\n\
The indices are themselves tuples (i,j) of row\n\
and column values.";

static PyObject *
LLMat_items(LLMatObject *a, PyObject *args) {
  PyObject *list;           /* the list that will hold the values */
  int i, j, k;
  int pos = 0;                /* position in list */
  double val;
    
  if (!PyArg_ParseTuple(args, "")) return NULL;
    
  if ((list = PyList_New(a->nnz)) == NULL) return NULL;
        
  for (i = 0; i < a->dim[0]; i ++) {
    for( k=a->root[i]; k != -1; k=a->link[k] ) {
      j = a->col[k];
      val = a->val[k];
      PyList_SET_ITEM(list, pos++, Py_BuildValue("((ii)d)", i, j, val));
    }
  }
  return list;
}

static char scale_doc[] = "A.scale(sigma)\n\
\n\
Scale each element in the matrix by the constant sigma.\n";

static PyObject *
LLMat_scale(LLMatObject *self, PyObject *args) {
  double sigma;
  int i, k;
  
  if (!PyArg_ParseTuple(args, "d", &sigma))
      return NULL;
  
  for (i = 0; i < self->dim[0]; i++)
    for( k=self->root[i]; k != -1; k=self->link[k] )
      self->val[k] *= sigma;
  
  Py_INCREF(Py_None); 
  return Py_None;
}


static char update_add_mask_doc[] = \
"A.update_add_mask(b, ind0, ind1, mask0, mask1)\n\
\n\
Update of global FEM matrix. Equivalent to:\n\
\n\
for i in range(len(ind0)):\n\
    for j in range(len(ind1)):\n\
        if mask0[i] and mask1[j]:\n\
            a[ind0[i],ind1[j]] += b[i,j]";

static PyObject *
LLMat_update_add_mask(LLMatObject *self, PyObject *args) {
  PyObject *bIn, *ind0In, *ind1In, *mask0In, *mask1In; 
  PyArrayObject *b, *ind0, *ind1, *mask0, *mask1;
  double v;
  int len0, len1, i, j, i1, j1, ldb;

  if (self->issym) {
    PyErr_SetString(SpMatrix_ErrorObject,
                    "Method not allowed for symmetric matrices");
    return NULL;
  }

  if (!PyArg_ParseTuple(args, "OOOOO",
                        &bIn, &ind0In, &ind1In, &mask0In, &mask1In)) 
    return NULL;

  b = (PyArrayObject *)PyArray_ContiguousFromObject(bIn, PyArray_DOUBLE, 2, 2);
  ind0 = (PyArrayObject *)PyArray_ContiguousFromObject(ind0In,PyArray_LONG,1,1);
  ind1 = (PyArrayObject *)PyArray_ContiguousFromObject(ind1In,PyArray_LONG,1,1);
  mask0 = (PyArrayObject *)PyArray_ContiguousFromObject(mask0In, PyArray_LONG,
                                                        1, 1);
  mask1 = (PyArrayObject *)PyArray_ContiguousFromObject(mask1In, PyArray_LONG,
                                                        1, 1);

  if (b==NULL || ind0==NULL || ind1==NULL || mask0==NULL || mask1==NULL)
    goto fail;

  len0 = ind0->dimensions[0];
  len1 = ind1->dimensions[0];

  /* validate array shapes */
  if (mask0->dimensions[0] != len0 || mask1->dimensions[0] != len1) {
    PyErr_SetString(PyExc_ValueError,
                    "Shapes of index and mask array do not match");
    goto fail;
  }
  if (b->dimensions[0] != len0 || b->dimensions[1] != len1) {
    PyErr_SetString(PyExc_ValueError,
                    "Shapes of input matrix and index arrays do not match");
    goto fail;
  }
  
  /* perform update add operation */
  ldb = b->dimensions[0];
  for (i = 0; i < len0; i ++) {
    if (((long *)mask0->data)[i]) {

      i1 = ((long *)ind0->data)[i];
      if (i1 < 0)
	i1 += self->dim[0];
      if (i1 < 0 || i1 >= self->dim[0]) {
	PyErr_SetString(PyExc_IndexError, "element of arg 2 out of range");
	goto fail;
      }
   
      for (j = 0; j < len1; j ++) {
	if (((long *)mask1->data)[j]) {

	  j1 = ((long *)ind1->data)[j];
	  if (j1 < 0)
	    j1 += self->dim[1];
	  if (j1 < 0 || j1 >= self->dim[1]) {
	    PyErr_SetString(PyExc_IndexError, "element of arg 3 out of range");
	    goto fail;
	  }

	  v = ((double *)b->data)[i + ldb*j];
	  if (SpMatrix_LLMatUpdateItemAdd(self, i1, j1, v) == -1)
	    goto fail;
	}
      }
    }
  }
  
  Py_DECREF(b);
  Py_DECREF(ind0);
  Py_DECREF(ind1);
  Py_DECREF(mask0);
  Py_DECREF(mask1);
  Py_INCREF(Py_None); 
  return Py_None;

 fail:
  Py_XDECREF(b);
  Py_XDECREF(ind0);
  Py_XDECREF(ind1);
  Py_XDECREF(mask0);
  Py_XDECREF(mask1);
  return NULL;
}

static char update_add_mask_sym_doc[] = "A.update_add_mask(b, ind, mask)\n\
\n\
Symmetric update of global FEM matrix. Equivalent to:\n\
\n\
for i in range(len(ind)):\n\
    for j in range(len(ind)):\n\
        if mask[i] and mask[i]:\n\
            i1 = ind[i]; j1 = ind[j]\n\
            if i >= j:\n\
                a[i1,j1] += b[i,j]\n\
\n\
Only operates on the lower triangle of A. Used for symmetric matrices.";

static PyObject *
LLMat_update_add_mask_sym(LLMatObject *self, PyObject *args) {
  PyObject *bIn, *indIn, *maskIn; 
  PyArrayObject *b, *ind, *mask;
  double v;
  int len, i, j, i1, j1, ldb;

  if (!PyArg_ParseTuple(args, "OOO", &bIn, &indIn, &maskIn)) return NULL;

  b = (PyArrayObject *)PyArray_ContiguousFromObject(bIn, PyArray_DOUBLE, 2, 2);
  ind = (PyArrayObject *)PyArray_ContiguousFromObject(indIn, PyArray_LONG, 1,1);
  mask = (PyArrayObject *)PyArray_ContiguousFromObject(maskIn,PyArray_LONG,1,1);

  if (b == NULL || ind == NULL || mask == NULL) goto fail;

  len = ind->dimensions[0];

  /* validate array shapes */
  if (mask->dimensions[0] != len) {
    PyErr_SetString(PyExc_ValueError,
                    "Shapes of index and mask array do not match");
    goto fail;
  }
  if (b->dimensions[0] != len || b->dimensions[1] != len) {
    PyErr_SetString(PyExc_ValueError,
                    "Shapes of input matrix and index arrays do not match");
    goto fail;
  }
  
  /* perform update add operation */
  ldb = b->dimensions[0];
  for (i = 0; i < len; i ++) {
    if (((long *)mask->data)[i]) {

      i1 = ((long *)ind->data)[i];
      if (i1 < 0)
	i1 += self->dim[0];
      if (i1 < 0 || i1 >= self->dim[0]) {
	PyErr_SetString(PyExc_IndexError, "element of arg 2 out of range");
	goto fail;
      }
   
      for (j = 0; j <= i; j ++) {
	if (((long *)mask->data)[j]) {

	  j1 = ((long *)ind->data)[j]; /* index check not necessary here */
	  if (j1 < 0)
	    j1 += self->dim[1];
	  v = ((double *)b->data)[i + ldb*j];

	  if (self->issym) {
	    /* symmetric matrix: update entry in lower triangle */
	    if (i1 > j1) {
	      if (SpMatrix_LLMatUpdateItemAdd(self, i1, j1, v) == -1)
		goto fail;
	    } else {
	      if (SpMatrix_LLMatUpdateItemAdd(self, j1, i1, v) == -1)
		goto fail;
	    }
	  } else {
	    /* non-symmetric matrix: update two entries if not on diagonal */
	    if (SpMatrix_LLMatUpdateItemAdd(self, i1, j1, v) == -1)
	      goto fail;
	    if (i1 != j1) {
	      if (SpMatrix_LLMatUpdateItemAdd(self, j1, i1, v) == -1)
		goto fail;
	    }
	  }

	}
      }
    }
  }
  
  Py_DECREF(b);
  Py_DECREF(ind);
  Py_DECREF(mask);
  Py_INCREF(Py_None); 
  return Py_None;

 fail:
  Py_XDECREF(b);
  Py_XDECREF(ind);
  Py_XDECREF(mask);
  return NULL;
}

static char LLMat_take_doc[] = "A.take(b,id1,id2)\n\
\n\
for i in range(len(b)):\n\
    b[i] = A[id1[i],id2[i]]";

static PyObject *
LLMat_take(LLMatObject *self, PyObject *args) {
  PyObject *bIn;
  PyObject *id1in = NULL;
  PyObject *id2in = NULL;
  PyArrayObject *b;
  PyArrayObject *id1 = NULL;
  PyArrayObject *id2 = NULL;
  int lenb,i;

  if (!PyArg_ParseTuple(args, "O|OO", &bIn,&id1in,&id2in)) return NULL;

  b = (PyArrayObject *)PyArray_ContiguousFromObject(bIn, PyArray_DOUBLE, 1, 1);
  if (b == NULL)
    goto fail;

  lenb = b->dimensions[0];

  if (id1in) {
    id1 = (PyArrayObject *)PyArray_ContiguousFromObject(id1in,PyArray_LONG,1,1);
    if (id1 == NULL) goto fail;
  }
  
  if (id2in) {
    id2 = (PyArrayObject *)PyArray_ContiguousFromObject(id2in,PyArray_LONG,1,1);
    if (id2 == NULL) goto fail;
  }
  
  if (lenb < 0 ) {
    PyErr_SetString(PyExc_IndexError, "vector b has a negative size");
    goto fail;
  }
  
  if (id1 && id1->dimensions[0] != lenb ) {
    PyErr_SetString(PyExc_IndexError, "id1 does not have the same size as b");
    goto fail;
  }

  if (id2 && id2->dimensions[0] != lenb ) {
    PyErr_SetString(PyExc_IndexError, "id2 does not have the same size as b");
    goto fail;
  }
  
  /*
  if (id1 != id2 && self->issym) {
    PyErr_SetString(SpMatrix_ErrorObject,
                    "Symmetric matrices require identical sets of indices");
    goto fail;
  }
  */
    
  /* Perform take operation */
  for( i = 0; i < lenb; i++ ) {
    int	i1, j1;
    if( id1 )
      i1 = ((long *) id1->data)[i];
    else
      i1 = i;

    if( id2 )
      j1 = ((long *) id2->data)[i];
    else
      j1 = i1;

    if( i1 > j1 || !self->issym ) /* Get entries as given */
      ((double *)b->data)[i] = SpMatrix_LLMatGetItem(self, i1, j1);
    else    /* Symmetric matrix: get entries from lower triangle */
      ((double *)b->data)[i] = SpMatrix_LLMatGetItem(self, j1, i1);
  }
      
  Py_DECREF(b);
  if (id1) {
    Py_DECREF(id1);
  }
  if (id2) {
    Py_DECREF(id2);
  }
  Py_INCREF(Py_None); 
  return Py_None;
  
 fail:
  Py_XDECREF(b);
  if (id1) {
    Py_XDECREF(id1);
  }
  if (id2) {
    Py_XDECREF(id2);
  }
  return NULL;
}

static char LLMat_put_doc[] = "a.put(b,id1,id2)\n\
\n\
for i in range(len(b)):\n\
    a[id1[i],id2[i]] = b[i]\n\n\
If b is a scalar, it has the same effect as the list [b, b, ..., b]\n\
If id1 is omitted, it is considered to be [1, 2, 3, ...].\n\
If id2 is omitted, it is considered equal to id1.\n";

static PyObject *
LLMat_put(LLMatObject *self, PyObject *args) {

  // For simplicity, we use iterators to parse Numpy arrays. We could gain in
  // efficiency by checking whether the array is contiguous or not. If it is,
  // a simple loop over the elements will be faster than the iterator.

  PyObject *bIn, *ind;
  PyObject *id1in = NULL;
  PyObject *id2in = NULL;
  //PyArrayObject *b = NULL;
  //PyArrayObject *id1 = NULL;
  //PyArrayObject *id2 = NULL;
  long lenb = 0, lenid1 = 0, lenid2 = 0, i;
  char b_is_scalar = 0, b_is_list = 0, id1_is_list = 0, id2_is_list = 0;
  double bval = 0.0;

  PyObject *iterator0 = NULL, *iterator1 = NULL, *iterator2 = NULL;

  if (!PyArg_ParseTuple(args, "O|OO", &bIn,&id1in,&id2in))
    return NULL;

  // Determine nature of value array b
  if( PyInt_Check(bIn) ) {                // b is an integer

    //printf("put: b is an Int\n");
    bval = (double)PyInt_AsLong(bIn);
    if( PyErr_Occurred() ) {
      PyErr_SetString(PyExc_TypeError,
                      "Could not convert int to double");
      goto fail;
    }
    b_is_scalar = 1;
    lenb = 1;

  } else if( PyLong_Check(bIn) ) {        // b in a long int

    bval = PyLong_AsDouble(bIn);
    b_is_scalar = 1;
    lenb = 1;

  } else if( PyFloat_Check(bIn) ) {       // b is a float

    //printf("put: b is a Float\n");
    bval = PyFloat_AsDouble(bIn);
    b_is_scalar = 1;
    lenb = 1;

  } else if( PyList_Check(bIn) ) {        // b is a list

    //printf("put: b is a list\n");
    lenb = (long)PyList_Size(bIn);
    b_is_list = 1;
    
  } else if( PyArray_Check(bIn) ) {       // b is an array

    //printf("put: b is an array\n");
    if( !PyArray_ISINTEGER(bIn) && !PyArray_ISFLOAT(bIn) ) {
      PyErr_SetString(PyExc_TypeError,
                      "Value array must be Int, Long or Float");
      goto fail;
    }
    iterator0 = PyArray_IterNew(bIn);
    lenb = (long)PyArray_DIM(bIn, 0);
    PyArray_ITER_RESET(iterator0);

  } else {

    PyErr_SetString(PyExc_TypeError,
                    "Values must be Int, Long, Float, list or Numpy array");
    goto fail;
  }

  if (lenb < 0 ) {
    PyErr_SetString(PyExc_IndexError, "vector b has a negative size");
    goto fail;
  }

  // Determine nature of first index, if given
  if (id1in) {

    if( PyList_Check(id1in) ) {           // id1 is a list

      //printf("put: id1in is a list\n");
      lenid1 = (long)PyList_Size(id1in);
      id1_is_list = 1;

    } else if( PyArray_Check(id1in) ) {   // id1 is a Numpy array

      //printf("put: id1in is an array\n");
      iterator1 = PyArray_IterNew(id1in); // id1 may not be contiguous
      lenid1 = (long)PyArray_DIM(id1in, 0);
      PyArray_ITER_RESET(iterator1);

    } else {

      PyErr_SetString(PyExc_TypeError,
                      "First index must be list or Numpy array");
      goto fail;
    }

    if( !b_is_scalar )
      if( lenid1 != lenb ) {
        PyErr_SetString(PyExc_IndexError,
                        "Not as many row indices as values");
        goto fail;
      }

    if( b_is_scalar ) lenb = lenid1;
  }
    
  // Determine nature of second index, if given
  if( id2in ) {

    if( PyList_Check(id2in) ) {           // id2 is a list

      //printf("put: id2in is a list\n");
      lenid2 = (long)PyList_Size(id2in);
      id2_is_list = 1;

    } else if( PyArray_Check(id2in) ) {   // id2 is a Numpy array

      //printf("put: id2in is an array\n");
      iterator2 = PyArray_IterNew(id2in); // id2 may not be contiguous
      lenid2 = (long)PyArray_DIM(id2in, 0);
      PyArray_ITER_RESET(iterator2);

    } else {

      PyErr_SetString(PyExc_TypeError,
                        "Second index must be list or Numpy array");
      goto fail;
    }

    if( !b_is_scalar ) {
      if( lenid2 != lenb ) {
        PyErr_SetString(PyExc_IndexError,
                        "Not as many column indices as values");
        goto fail;
      }
    }

    if( id1in ) {
      if( lenid1 != lenid2 ) {
        PyErr_SetString(PyExc_IndexError,
                        "Not as many row indices as column indices");
        goto fail;
      }
    }

    if( b_is_scalar ) lenb = lenid2;
  }

  // Perform put operation
  for( i = 0; i < lenb; i++ ) {
    long i1, j1;
    
    i1 = i;
    if( id1in ) {
      if( id1_is_list ) {
        ind = PyList_GetItem(id1in, (Py_ssize_t)i);
        if( PyInt_Check(ind) )
          i1 = PyInt_AS_LONG(ind);
        else {
          PyErr_SetString(PyExc_ValueError, "Invalid list item");
          return NULL;
        }
      } else {
        i1 = *(long*)PyArray_ITER_DATA(iterator1);
        PyArray_ITER_NEXT(iterator1);
      }
    }

    j1 = i1;
    if( id2in ) {
      if( id2_is_list ) {
        ind = PyList_GetItem(id2in, (Py_ssize_t)i);
        if( PyInt_Check(ind) )
          j1 = PyInt_AS_LONG(ind);
        else {
          PyErr_SetString(PyExc_ValueError, "Invalid list item");
          return NULL;
        }
      } else {
        j1 = *(long*)PyArray_ITER_DATA(iterator2);
        PyArray_ITER_NEXT(iterator2);
      }
    }

    if( !b_is_scalar ) {
      if( b_is_list ) {
        ind = PyList_GetItem(bIn, (Py_ssize_t)i);
        if( PyInt_Check(ind) )
          bval = (double)PyInt_AS_LONG(ind);
        else if( PyFloat_Check(ind) )
          bval = PyFloat_AsDouble(ind);
        else {
          PyErr_SetString(PyExc_ValueError, "Invalid list item");
          return NULL;
        }
      } else {
        // Convert value to double appropriately
        if( PyArray_ISINTEGER(bIn) )
          bval = (double)(*(long*)(PyArray_ITER_DATA(iterator0)));
        else  // float
          bval = *(double*)(PyArray_ITER_DATA(iterator0));
        PyArray_ITER_NEXT(iterator0);
      }
    }

    //printf(" %g  --> (%ld,%ld)\n", bval, i1, j1);

    if (i1 > j1 || !self->issym) { /* Update entries as given */
      //printf("mat[%d,%d] <- %g\n", i1, j1, bval);
      if( SpMatrix_LLMatSetItem(self, i1, j1, bval) == -1 )
        goto fail;
    } else { /* Symmetric matrix: update entries in lower triangle */
      if( SpMatrix_LLMatSetItem(self, j1, i1, bval) == -1 )
        goto fail;
    }
  }
  
  /*
    if( !b_is_scalar ) {
    Py_DECREF(b);
    }
    if (id1) {
    Py_DECREF(id1);
    }
    if (id2) {
    Py_DECREF(id2);
    }
  */
  Py_XDECREF(iterator0);
  Py_XDECREF(iterator1);
  Py_XDECREF(iterator2);
  Py_INCREF(Py_None); 
  return Py_None;
  
 fail:
  /*
    if( !b_is_scalar ) {
    Py_XDECREF(b);
    }
    if (id1) {
    Py_XDECREF(id1);
    }
    if (id2) {
    Py_XDECREF(id2);
    }
  */
  Py_XDECREF(iterator0);
  Py_XDECREF(iterator1);
  Py_XDECREF(iterator2);
  return NULL;
}

static char LLMat_delete_rows_doc[] = 
"Delete rows from matrix (inplace). The rows to be deleted are specified by the mask array.\n\
\n\
Arguments:\n\
\n\
  mask: A 1D integer NumPy array. If mask[i] == 0, then row i is deleted,\n\
        otherwise row i is kept.\n\
\n\
This method may not be applied to a matrix in symmetric format.";

static PyObject*
LLMat_delete_rows(LLMatObject *self, PyObject* args){
  PyArrayObject *maskObj;
  int newm, newnnz;
  int act, row;
  
  if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &maskObj)) 
    return NULL;
  if (maskObj->nd != 1 || maskObj->descr->type_num != PyArray_LONG || maskObj->dimensions[0] != self->dim[0]) {
    PyErr_SetString(PyExc_ValueError, "mask must be a 1D integer NumPy array of appropriate length");
    return NULL;
  }
  if (self->issym) {
    PyErr_SetString(SpMatrix_ErrorObject, "method not allowed for symmetric matrices");
    return NULL;
  }

  /* Delete the rows to be cancelled by rearranging the row */
  /* array. After having done so, newdim is the new matrix dim. */
  newm = 0;
  newnnz = self->nnz;
  for(row = 0; row < self->dim[0]; row ++){
    if (*(int *)(maskObj->data + row*maskObj->strides[0]) != 0){ /* This row has to be kept */
      self->root[newm] = self->root[row];
      newm ++;
    } else {			/* row let out; update free list */
      act = self->root[row];
      if(act != -1){		/* only do smth. for non-empty rows */
	newnnz --;
	while(self->link[act] != -1) { /* Walk to the end of the list */
	  act = self->link[act];
	  newnnz --;
	}
	self->link[act] = self->free;	/* Attach end of row to free list */
	self->free = self->root[row];	/* Start free list where row began */
      }
    }
  }

  /* Set the new values */
  self->dim[0] = newm;
  self->nnz = newnnz;

  Py_INCREF(Py_None); 
  return Py_None;
}

static char LLMat_delete_cols_doc[] = 
"Delete columns from matrix (inplace). The columns to be deleted are\n\
specified by the mask array.\n\
\n\
Arguments:\n\
\n\
  mask: A 1D integer NumPy array. If mask[i] == 0, then column i is deleted,\n\
        otherwise column i is kept.";

static PyObject*
LLMat_delete_cols(LLMatObject *self, PyObject* args){
  PyArrayObject *maskObj;
  int newn, newnnz;
  int act, old, col, row;
  int *shift;
  
  if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &maskObj)) 
    return NULL;
  if (maskObj->nd != 1 || maskObj->descr->type_num != PyArray_LONG || maskObj->dimensions[0] != self->dim[1]) {
    PyErr_SetString(PyExc_ValueError,
                    "mask must be a 1D integer NumPy array of appropriate length");
    return NULL;
  }
  if (self->issym) {
    PyErr_SetString(SpMatrix_ErrorObject,
                    "method not allowed for symmetric matrices");
    return NULL;
  }

#define MASK(i) *(long *)(maskObj->data + (i)*maskObj->strides[0])

  /* Allocate column shift vector (after deletion col[i] is at */
  /* col[i] - shift[i]). */
  shift = (int*)malloc((self->dim[1])*sizeof(int));
  newn = self->dim[1];
  if (MASK(0)) shift[0] = 0; else {shift[0] = 1; newn --;}
  for (col = 1; col < self->dim[1]; col++){
    if (MASK(col)) 
      shift[col] = shift[col-1]; 
    else 
      {shift[col] = shift[col-1]+1; newn --; }
  }
    
  /* Deleteting columns in the remainig rows */
  newnnz = self->nnz;
  for(row = 0; row < self->dim[0]; row ++) {
    old = -1; act = self->root[row];
    while (act != -1){
      if (MASK(self->col[act])) {	      // Keep this column
	self->col[act] -= shift[self->col[act]];
	old = act; act = self->link[act];
      } else {				      // Drop the column
	newnnz--;                              
	if (self->root[row] == act) {	      // Special case: first row element
	  self->root[row] = self->link[act];
	  old = act; act = self->link[act];
	  self->link[old] = self->free;	      // Append element into freelist
	  self->free = old;
	} else {			      // Regular case: element inbetween
	  act = self->link[act];
	  self->link[self->link[old]] = self->free;
	  self->free = self->link[old];
	  self->link[old] = act;	      // Append element into freelist
	}
      }
    }
  }

  /* Set the new values */
  self->dim[1] = newn;
  self->nnz = newnnz;

  /* clean up */
  free(shift);

  Py_INCREF(Py_None); 
  return Py_None;

#undef MASK
}

static char LLMat_delete_rowcols_doc[] = 
"Delete rows and columns from matrix (inplace). The rows and columns to be deleted are\n\
specified by the mask array.\n\
\n\
Arguments:\n\
\n\
  mask: A 1D integer NumPy array. If mask[i] == 0, then row and column i are deleted,\n\
        otherwise row and column i are kept.";

static PyObject*
LLMat_delete_rowcols(LLMatObject *self, PyObject* args){
  PyArrayObject *maskObj;
  int newn, newm, newnnz;
  int act, old, col, row;
  int *shift;
  
  if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &maskObj)) 
    return NULL;
  if (maskObj->nd != 1 || maskObj->descr->type_num != PyArray_LONG || maskObj->dimensions[0] != self->dim[0]) {
    PyErr_SetString(PyExc_ValueError, "mask must be a 1D integer NumPy array of appropriate length");
    return NULL;
  }
  if (self->dim[0] != self->dim[1]) {
    PyErr_SetString(SpMatrix_ErrorObject, "method only allowed for square matrices");
    return NULL;
  }

#define MASK(i) *(long *)(maskObj->data + (i)*maskObj->strides[0])

  /* Delete the rows to be cancelled by rearranging the row */
  /* array. After having done so, newdim is the new matrix dim. */
  newm = 0;
  newnnz = self->nnz;
  for(row = 0; row < self->dim[0]; row ++){
    if (MASK(row)){			      // This row has to be kept
      self->root[newm] = self->root[row];
      newm ++;
    } else {				      // row let out; update free list
      act = self->root[row];
      if(act != -1){			      // only do sth for non-empty rows
	newnnz --;
	while(self->link[act] != -1) {	      // Walk to the end of the list
	  act = self->link[act];
	  newnnz --;
	}
	self->link[act] = self->free;	      // Attach end of row to free list
	self->free = self->root[row];	      // Start free list where row began
      }
    }
  }

  /* Set the new values */
  self->dim[0] = newm;
  self->nnz = newnnz;

  /* Allocate column shift vector (after deletion col[i] is at */
  /* col[i] - shift[i]). */
  shift = (int*)malloc((self->dim[1])*sizeof(int));
  newn = self->dim[1];
  if (MASK(0)) shift[0] = 0; else {shift[0] = 1; newn --;}
  for (col = 1; col < self->dim[1]; col++){
    if (MASK(col)) 
      shift[col] = shift[col-1]; 
    else 
      {shift[col] = shift[col-1]+1; newn --; }
  }
    
  /* Deleteting columns in the remainig rows */
  newnnz = self->nnz;
  for(row = 0; row < self->dim[0]; row ++) {
    old = -1; act = self->root[row];
    while (act != -1){
      if (MASK(self->col[act])) {	       /* Keep this column */
	self->col[act] -= shift[self->col[act]];
	old = act; act = self->link[act];
      } else {				       /* Drop the column */
	newnnz--;                              
	if (self->root[row] == act) {	      // Special case: first row element
	  self->root[row] = self->link[act];
	  old = act; act = self->link[act];
	  self->link[old] = self->free;	      // Append element into freelist
	  self->free = old;
	} else {			      // Regular case: element inbetween
	  act = self->link[act];
	  self->link[self->link[old]] = self->free;
	  self->free = self->link[old];
	  self->link[old] = act;	      // Append element into freelist
	}
      }
    }
  }

  /* Set the new values */
  self->dim[1] = newn;
  self->nnz = newnnz;

  /* clean up */
  free(shift);

  Py_INCREF(Py_None); 
  return Py_None;

#undef MASK
}

static char LLMat_Find_Doc[] = "Get LL matrix in coord format (val,irow,jcol).";

static PyObject *LLMat_Find( LLMatObject *self, PyObject *args ) {

  /* Convert an LL matrix into coordinate format */

  PyArrayObject *a_row, *a_col, *a_val; /* Matrix in coordinate format */
  npy_intp       dmat[1];               /* Dimension descriptor */
  int           *pi, *pj;             /* Intermediate pointers to matrix data */
  double        *pv;
  int            i, k, elem;

  dmat[0] = (npy_intp)(self->nnz);

  /* Allocate numarrays */
  a_row = (PyArrayObject *)PyArray_SimpleNew( 1, dmat, NPY_INT );
  a_col = (PyArrayObject *)PyArray_SimpleNew( 1, dmat, NPY_INT );
  a_val = (PyArrayObject *)PyArray_SimpleNew( 1, dmat, NPY_DOUBLE );

  pi = (int *)a_row->data;
  pj = (int *)a_col->data;
  pv = (double *)a_val->data;

  elem = 0;
  for( i = 0; i < self->dim[0]; i++ ) {
    for( k = self->root[i]; k != -1; k = self->link[k] ) {
      pi[ elem ] = i;
      pj[ elem ] = self->col[k];
      pv[ elem ] = self->val[k];
      elem++;
    }
  }

  return Py_BuildValue( "OOO",
                        PyArray_Return( a_val ),
                        PyArray_Return( a_row ),
                        PyArray_Return( a_col ) );
}

/*********************************/
/*  O b j e c t   m e t h o d s  */
/*********************************/

PyMethodDef LLMat_methods[] = {
  {"matvec",        (PyCFunction)LLMat_matvec,        METH_VARARGS, LLMat_matvec_doc},
  {"matvec_transp", (PyCFunction)LLMat_matvec_transp, METH_VARARGS, LLMat_matvec_transp_doc},
  {"to_csr",        (PyCFunction)LLMat_to_csr,        METH_VARARGS, to_csr_doc},
  {"to_sss",        (PyCFunction)LLMat_to_sss,        METH_VARARGS, to_sss_doc},
  {"generalize",    (PyCFunction)LLMat_generalize,    METH_VARARGS, LLMat_generalize_doc},
  {"compress",      (PyCFunction)LLMat_compress,      METH_VARARGS, LLMat_compress_doc},
  {"export_mtx",    (PyCFunction)LLMat_export_mtx,    METH_VARARGS, export_mtx_doc},
  {"copy",          (PyCFunction)LLMat_copy,          METH_VARARGS, copy_doc},
  {"norm",          (PyCFunction)LLMat_norm,          METH_VARARGS, LLMat_norm_doc},
  {"shift",         (PyCFunction)LLMat_shift,         METH_VARARGS, shift_doc},
  {"scale",         (PyCFunction)LLMat_scale,         METH_VARARGS, scale_doc},
  {"col_scale",     (PyCFunction)LLMat_col_scale,     METH_VARARGS, col_scale_doc},
  {"row_scale",     (PyCFunction)LLMat_row_scale,     METH_VARARGS, row_scale_doc},
  {"keys",          (PyCFunction)LLMat_keys,          METH_VARARGS, keys_doc},
  {"values",        (PyCFunction)LLMat_values,        METH_VARARGS, values_doc},
  {"items",         (PyCFunction)LLMat_items,         METH_VARARGS, items_doc},
  {"put",           (PyCFunction)LLMat_put,           METH_VARARGS, LLMat_put_doc},
  {"take",          (PyCFunction)LLMat_take,          METH_VARARGS, LLMat_take_doc},
  { "find",         (PyCFunction)LLMat_Find, METH_VARARGS, LLMat_Find_Doc  },
  {"update_add_mask", (PyCFunction)LLMat_update_add_mask, METH_VARARGS, update_add_mask_doc},
  {"update_add_mask_sym", (PyCFunction)LLMat_update_add_mask_sym, METH_VARARGS, update_add_mask_sym_doc},
  {"delete_rows",   (PyCFunction)LLMat_delete_rows,     METH_VARARGS, LLMat_delete_rows_doc},
  {"delete_cols",   (PyCFunction)LLMat_delete_cols,     METH_VARARGS, LLMat_delete_cols_doc},
  {"delete_rowcols", (PyCFunction)LLMat_delete_rowcols,  METH_VARARGS, LLMat_delete_rowcols_doc},
  {"update_add_at", (PyCFunction)LLMat_update_add_at,  METH_VARARGS, update_add_at_doc},

  {NULL, NULL}			/* sentinel */
};

/*****************************************/
/*  L L M a t   t y p e   m e t h o d s  */
/*****************************************/

static void
LLMatType_dealloc(LLMatObject *a)
{
  PyMem_DEL(a->root);
  PyMem_DEL(a->val);
  PyMem_DEL(a->col);
  PyMem_DEL(a->link);
  PyObject_Del(a);
}

static int
LLMatType_print(LLMatObject *a, FILE *fp, int flags)
{
  int i, k, first = 1;
  char *symStr;

  if (a->issym)
    symStr = "symmetric";
  else
    symStr = "general";

  if (a->dim[1] <= PPRINT_COL_THRESH && a->dim[0] <= PPRINT_ROW_THRESH) {
    double *mat;
    int j;
    double val;
    mat = (double *)malloc(a->dim[0]*a->dim[1] * sizeof(double));
    if (mat == NULL) {
      PyErr_NoMemory();
      return -1;
    }
    fprintf(fp, "ll_mat(%s, [%d,%d]):\n", symStr, a->dim[0], a->dim[1]);
    for (i = 0; i < a->dim[0]; i ++) {
      for (j = 0; j < a->dim[1]; j ++)
	mat[i*a->dim[1] + j] = 0.0;
      k = a->root[i];
      while (k != -1) {
	mat[(i*a->dim[1])+a->col[k]] = a->val[k];
	k = a->link[k];
      }
    }

    for (i = 0; i < a->dim[0]; i ++) {
      for (j = 0; j < a->dim[1]; j ++) {
	val = mat[(i*a->dim[1])+j];
	if (val != 0.0) {
	  int exp = (int)log10(fabs(val));
	  if (abs(exp) <= 4) {
	    if (exp < 0)
	      fprintf(fp, "%9.*f ", 6, val);
	    else
	      fprintf(fp, "%9.*f ", 6-exp, val);
	  } else
	    fprintf(fp, "%9.1e ", val);
	}
	else
	  if (!(a->issym) || i > j)
	    fprintf(fp, " -------- ");
      }
      fprintf(fp, "\n");
    }
    free(mat);

  } else {

  if (a->nnz == 0) {
    fprintf(fp, "ll_mat(%s, [%d,%d])", symStr, a->dim[0], a->dim[1]);
    return 0;
  }
  fprintf(fp, "ll_mat(%s, [%d,%d], [", symStr, a->dim[0], a->dim[1]);
  for (i = 0; i < a->dim[0]; i ++) {
    k = a->root[i];
    while (k != -1) {
      if (!first)
	fprintf(fp, ", ");
      first = 0;
      fprintf(fp, "(%d,%d): %g", i, a->col[k], a->val[k]);
      k = a->link[k];
    }
  }
  fprintf(fp, "])");
  }
  return 0;
}

static PyObject *
LLMatType_getattr(LLMatObject *self, char *name)
{
  if (strcmp(name, "shape") == 0)
    return Py_BuildValue("(i,i)", self->dim[0], self->dim[1]);
  if (strcmp(name, "nnz") == 0)
    return PyInt_FromLong(self->nnz);
  if (strcmp(name, "issym") == 0)
    return PyInt_FromLong(self->issym);
  if (strcmp(name, "__members__") == 0) {
    char *members[] = {"shape", "nnz", "issym"};
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
  return Py_FindMethod(LLMat_methods, (PyObject *)self, name);
}
/***********************************************************************
 * mapping functions
 */

/** LLMat_length - number of items in mapping
 *    == number of matrix entries
 */
static int LLMat_length(LLMatObject *self) {
  return self->dim[0] * self->dim[1];
}

/** LLMat_subscript
 *    Called when treating array object like a mapping. This is used to
 *    implement two-dimensional idices, e.g. A[i,j] or A[i1:i2:i3,j1:j2]
 */

static PyObject *LLMat_subscript(LLMatObject *self, PyObject *index) {

  PyObject *index0, *index1;

  // Check that input is a double index
  if( !PySequence_Check(index) ) {
    PyErr_SetString(PyExc_IndexError, "Index must be a sequence");
    Py_INCREF(Py_None);
    return Py_None;
  }
  if( PySequence_Length(index) != 2 ) {
    PyErr_SetString(PyExc_IndexError, "There must be exactly two indices");
    Py_INCREF(Py_None);
    return Py_None;
  }

  // Parse first index
  if( !(index0 = PySequence_GetItem(index, 0)) ) {
    PyErr_SetString(PyExc_IndexError, "First index is invalid");
    Py_INCREF(Py_None);
    return Py_None;
  }

  // Parse second index
  if( !(index1 = PySequence_GetItem(index, 1)) ) {
    PyErr_SetString(PyExc_IndexError, "Second index is invalid");
    Py_INCREF(Py_None);
    return Py_None;
  }

  // Return submatrix
  return getSubMatrix_FromList(self, index0, index1);

}

/** LLMat_ass_subscript
 *    Called when treating array object like a mapping. This is used
 *    implement two-dimensional indices, e.g. A[::2,1:5]
 */

static int
LLMat_ass_subscript(LLMatObject *self, PyObject *index, PyObject *value ) {

  PyObject *index0, *index1;

  // Check that input is a double index
  if( !PySequence_Check(index) ) {
    PyErr_SetString(PyExc_IndexError, "Index must be a sequence");
    return -1;
  }
  if( PySequence_Length(index) != 2 ) {
    PyErr_SetString(PyExc_IndexError, "There must be exactly two indices");
    return -1;
  }

  // Parse first index
  if( !(index0 = PySequence_GetItem(index, 0)) ) {
    PyErr_SetString(PyExc_IndexError, "First index is invalid");
    return -1;
  }

  // Parse second index
  if( !(index1 = PySequence_GetItem(index, 1)) ) {
    PyErr_SetString(PyExc_IndexError, "Second index is invalid");
    return -1;
  }

  // Assign a submatrix
  //return setSubMatrix_FromList(self, value, index0, index1);
  setSubMatrix_FromList(self, value, index0, index1);
  if( PyErr_Occurred() )
    return -1;
  return 0;
}

static PyMappingMethods LLMat_as_mapping = {
#ifdef LENFUNC_OK
  (lenfunc)LLMat_length,              /*mp_length*/
#else
  (inquiry)LLMat_length,	      /*mp_length*/
#endif
  (binaryfunc)LLMat_subscript,        /*mp_subscript*/
  (objobjargproc)LLMat_ass_subscript, /*mp_ass_subscript*/
};

/*************************************/
/*  L L M a t T y p e   o b j e c t  */
/*************************************/

static PyTypeObject LLMatType = {
  PyObject_HEAD_INIT(NULL)
  0,                              /* ob_size */
  "ll_mat",                       /* tp_name */
  sizeof(LLMatObject),            /* tp_basicsize */
  0,                              /* tp_itemsize */
  (destructor)LLMatType_dealloc,  /* tp_dealloc */
  (printfunc)LLMatType_print,	  /* tp_print */
  (getattrfunc)LLMatType_getattr, /* tp_getattr */
  0,				  /* tp_setattr */
  0,				  /* tp_compare */
  0,				  /* tp_repr */
  0,				  /* tp_as_number */
  0,				  /* tp_as_sequence */
  &LLMat_as_mapping,		  /* tp_as_mapping */
  0,				  /* tp_hash */
  0,                              /* tp_call */
  0,                              /* tp_str */
  0,                              /* tp_getattro */
  0,                              /* tp_setattro */
  0,                              /* tp_as_buffer */
  Py_TPFLAGS_DEFAULT,             /* tp_flags */
  "LLMat objects",                /* tp_doc */
};

/*************************************************/
/*  M i s c .   h e l p e r   f u n c t i o n s  */
/*************************************************/

static PyObject *
SpMatrix_NewLLMatObject(int dim[], int sym, int sizeHint) {
  int i;
  LLMatObject *op;

  if (dim[0] < 0 || dim[1] < 0) {
    PyErr_SetString(PyExc_ValueError, "matrix dimension must be non-negative");
    return NULL;
  }
  if (sizeHint < 1)
    sizeHint = 1;

  /* create new SparseArrayt object */
  op = PyObject_New(LLMatObject, &LLMatType);
  if (op == NULL)
    return PyErr_NoMemory();

  op->root = NULL;
  op->val = NULL;
  op->col = NULL;
  op->link = NULL;

  /* allocate ob_val and on_idx arrays */
  op->root = PyMem_New(int, dim[0]);
  if (op->root == NULL)
    goto fail;
  op->val = PyMem_New(double, sizeHint);
  if (op->val == NULL)
    goto fail;
  op->col = PyMem_New(int, sizeHint);
  if (op->col == NULL)
    goto fail;
  op->link = PyMem_New(int, sizeHint);
  if (op->link == NULL)
    goto fail;

  /* initialize rest of fields */
  for (i = 0; i < dim[0]; i ++)
    op->root[i] = -1;
  op->dim[0] = dim[0];
  op->dim[1] = dim[1];
  op->issym = sym;
  op->nnz = 0;
  op->nalloc = sizeHint;
  op->free = -1;

  return (PyObject *) op;

 fail:
    PyMem_Del(op->link);    
    PyMem_Del(op->col);    
    PyMem_Del(op->val);    
    PyMem_Del(op->root);    
    PyObject_Del(op);
    return PyErr_NoMemory();
}

static PyObject*
LLMat_from_mtx(PyObject *module, PyObject *args) {
  LLMatObject *self = NULL;
  char *fileName;
  MM_typecode matcode;
  int dim[2], nz;
  FILE *f;
  int ret, i;
  double val;
  int row, col;

  if (!PyArg_ParseTuple(args, "s", &fileName))
    return NULL;
  
  /* open file */
  f = fopen(fileName, "r");
  if (f == NULL)
    return PyErr_SetFromErrno(PyExc_IOError);

  /* read MTX header */
  ret = mm_read_banner(f, matcode);
  if (ret != 0) {
    PyErr_SetString(PyExc_IOError, "error reading MTX file header");
    goto fail;
  }
  if (!(mm_is_real(matcode) && mm_is_matrix(matcode) &&
	mm_is_sparse(matcode))) {
    PyErr_SetString(SpMatrix_ErrorObject, "must be real, sparse matrix");
    goto fail;
  }
  ret = mm_read_mtx_crd_size(f, dim, dim+1, &nz);
  if (ret != 0) {
    PyErr_SetString(PyExc_IOError, "error reading MTX file size information");
    goto fail;
  }

  /* allocate matrix object */
  self = (LLMatObject *)SpMatrix_NewLLMatObject(dim, mm_is_symmetric(matcode), nz);
  if (self == NULL)
    goto fail;

  for (i = 0; i < nz; i ++) {
    ret = fscanf(f, "%d %d %lg\n", &row, &col, &val);
    if (ret != 3) {
      PyErr_SetString(PyExc_IOError, "error reading MTX file data");
      goto fail;
    }
    row --; col --;
    if (!(0 <= row && row < dim[0] && 0 <= col && col < dim[1])) {
      PyErr_SetString(PyExc_IndexError, //SpMatrix_ErrorObject,
                      "matrix indices out of range");
      fclose(f);
      return NULL;
    }
    ret = SpMatrix_LLMatSetItem(self, row, col, val);
    if (ret)
      goto fail;
  }
  fclose(f);
  return (PyObject *)self;

 fail:
  fclose(f);
  Py_XDECREF(self);
  return NULL;
}

char LLMat_matrixmultiply_doc[] = "matrixmultiply(A, B)\n\
\n\
Returns a new ll_mat object representing the matrix A*B";

static PyObject *
LLMat_matrixmultiply(PyObject *self, PyObject *args)
{
  int sizeHint = 1000;
  LLMatObject *matA, *matB, *matC;
  int dimC[2];
  int symCode, ret;
  
  if (!PyArg_ParseTuple(args, "O!O!", &LLMatType, &matA, &LLMatType, &matB))
    return NULL;

  /* matrix dimensions
   */
  dimC[0] = matA->dim[0];
  dimC[1] = matB->dim[1];

  if (matA->dim[1] != matB->dim[0]) {
    PyErr_SetString(PyExc_ValueError, "matrix dimensions must agree");
    return NULL;
  }

  /* create result object
   */
  matC = (LLMatObject *)SpMatrix_NewLLMatObject(dimC, 0, sizeHint);
  if (matC == NULL)
    return NULL;

  symCode = matB->issym << 1 | matA->issym;
  if (symCode == 0) {
    /* unsym-unsym multiplication
     */

#if !OPT_MATMATMUL
    double valA;
    int iA, jA, kA, kB; 

    for (iA = 0; iA < matA->dim[0]; iA ++) {
      kA = matA->root[iA];
      while (kA != -1) {
	valA = matA->val[kA];
	jA = matA->col[kA];
	kA = matA->link[kA];
	
	/* add jA-th row of B to iA-th row of C */
	kB = matB->root[jA];
	while (kB != -1) {
	  ret = SpMatrix_LLMatUpdateItemAdd(matC, iA, matB->col[kB], valA*matB->val[kB]);
	  if (ret == -1)
	    goto fail;
	  kB = matB->link[kB];
	}
      }

    }
#else
    int *tmpage = NULL;
    int *tmpind = NULL;
    int *tmpcol = NULL;
    double *tmpval = NULL;
    int tmpsize;
    int nxttmp;
    int row;
    int indA, colA, colB, dummy, indB;
    double valA;
    
    tmpsize = 5; nxttmp = -1;
    tmpage = (int*)malloc(matB->dim[1] * sizeof(int));
    tmpind = (int*)malloc(matB->dim[1] * sizeof(int));
    tmpcol = (int*)malloc(tmpsize * sizeof(int));
    tmpval = (double*)malloc(tmpsize * sizeof(double));
    if (tmpage == NULL || tmpind == NULL || tmpcol == NULL ||tmpval == NULL) {
      PyErr_NoMemory();
      goto fail_unsym_unsym;
    }

    /* main loop */
    
    for(row=0; row < matB->dim[1]; row++){ tmpage[row] = -1;}

    /* Go through each row of A and perform necessary computations */
    for(row=0; row < matA->dim[0]; row++) {
      indA = matA->root[row];         // Pick first entry of A[row,:]
      while(indA != -1){       // As long as there is an element in A[row,:] ...
	colA = matA->col[indA];       // ... get its column number ...
	valA = matA->val[indA];       // ... and value ...
	
	indB = matB->root[colA];       // colA is equivalent to rowB!
	while(indB != -1){
	  colB = matB->col[indB];
	  
	  if(tmpage[colB] != row){         // This column never appeared so far
	    nxttmp++;
	    tmpage[colB]  = row;   
	    tmpind[colB]  = nxttmp;
	    
	    if(nxttmp >= tmpsize){          // If tmp storage too small, realloc
	      tmpsize = (int)((tmpsize*12)/10)+1;
	      tmpcol = (int*)realloc(tmpcol, tmpsize * sizeof(int));
	      tmpval = (double*)realloc(tmpval, tmpsize * sizeof(double));
	      if (tmpcol == NULL ||tmpval == NULL) {
		PyErr_NoMemory();
		goto fail_unsym_unsym;
	      }
	    }
	    
	    tmpcol[nxttmp] = colB;
	    tmpval[nxttmp] = valA * matB->val[indB];
	  }else{                   // This column appeared at least once already
	    dummy = tmpind[colB];
	    tmpval[dummy] += valA * matB->val[indB];
	  }
	  
	  indB = matB->link[indB];
	}
	indA = matA->link[indA];
      }
	
      /* All the new values for rowC = rowA have now to be filled in */
      /* into the matrix C */
      for(dummy=0; dummy<=nxttmp; dummy++) {
	if (SpMatrix_LLMatSetItem(matC,row,tmpcol[dummy],tmpval[dummy]))
	  goto fail_unsym_unsym;
      }
      
      nxttmp=-1; /* For the next row of A we need a "fresh" tmp storage */
    }
    /* Get the memory back ... */
    free(tmpage);
    free(tmpind);
    free(tmpcol);
    free(tmpval);
    return (PyObject *)matC;

  fail_unsym_unsym:
    free(tmpage);
    free(tmpind);
    free(tmpcol);
    free(tmpval);
    goto fail;
#endif

  } else if (symCode == 1) {
    
    /* sym-unsym multiplication
     */
    double valA;
    int iA, jA, kA, kB;
    
    for (iA = 0; iA < matA->dim[0]; iA ++) {
      kA = matA->root[iA];
      while (kA != -1) {
	valA = matA->val[kA];
	jA = matA->col[kA];
	kA = matA->link[kA];

	/* add jA-th row of B to iA-th row of C */
	kB = matB->root[jA];
	while (kB != -1) {
	  ret = SpMatrix_LLMatUpdateItemAdd(matC, iA, matB->col[kB], valA*matB->val[kB]);
	  if (ret == -1)
	    goto fail;
	  kB = matB->link[kB];
	}
	
	if (iA == jA)
	  continue;
	
	/* add iA-th row of B to jA-th row of C */
	kB = matB->root[iA];
	while (kB != -1) {
	  ret = SpMatrix_LLMatUpdateItemAdd(matC, jA, matB->col[kB], valA*matB->val[kB]);
	  if (ret == -1)
	    goto fail;
	  kB = matB->link[kB];
	}	
      }
    }

  } else if (symCode == 2) {

    /* unsym-sym multiplication
     */

    PyErr_SetString(PyExc_NotImplementedError, "multiply of an unsymmetric and a symmetric matrix not supported");
    goto fail;

  } else {

    /* sym-sym multiplication
     */

    PyErr_SetString(PyExc_NotImplementedError, "multiply of two symmetric matrices not supported");
    goto fail;

  }
  
  return (PyObject *)matC;

 fail:
  Py_DECREF(matC);
  return NULL;
}

static char LLMat_dot_doc[] = "dot(A, B)\n\
\n\
Returns a new ll_mat object representing the matrix transpose(A)*B";

static PyObject *LLMat_dot(PyObject *self, PyObject *args) {

  int sizeHint = 1000;
  LLMatObject *matA, *matB, *matC;
  int dimC[2];
  double valA;
  int iA, kA, iC, kB, ret;
  
  if (!PyArg_ParseTuple(args, "O!O!", &LLMatType, &matA, &LLMatType, &matB))
    return NULL;

  dimC[0] = matA->dim[1];
  dimC[1] = matB->dim[1];

  if (matA->dim[0] != matB->dim[0]) {
    PyErr_SetString(PyExc_ValueError, "matrix dimensions must agree");
    return NULL;
  }

  if (matA->issym || matB->issym) {
    PyErr_SetString(PyExc_NotImplementedError,
                    "ddot operation with symmetric matrices not supported");
    return NULL;
  }

  matC = (LLMatObject *)SpMatrix_NewLLMatObject(dimC, 0, sizeHint);
  if (matC == NULL)
    return NULL;

  for( iA = 0; iA < matA->dim[0]; iA++ ) {
    for( kA = matA->root[iA]; kA != -1; kA = matA->link[kA] ) {
      valA = matA->val[kA];
      iC = matA->col[kA];
      for( kB = matB->root[iA]; kB != -1; kB = matB->link[kB] ) {
	ret = SpMatrix_LLMatUpdateItemAdd(matC, iC, matB->col[kB],
                                          valA*matB->val[kB]);
	if (ret == -1) goto fail;
      }
    }
  }
  return (PyObject *)matC;

 fail:
  Py_DECREF(matC);
  return NULL;
}

/* For backward compatibility. This is still called by sss_mat. */

static int 
LLMat_parse_index(PyObject *op, int dim[],
		  int *start0, int *stop0, int *step0, int *len0,
		  int *start1, int *stop1, int *step1, int *len1) {


  PyErr_SetString(PyExc_IndexError, "Not yet transitioned to fancy indexing for SSS matrices");
  return -1;
}
