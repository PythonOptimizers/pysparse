#include <float.h>

#include "Python.h"
#include "fortran.h"

#define PY_ARRAY_UNIQUE_SYMBOL sparslab
#include "Numeric/arrayobject.h"

#include "spmatrix.h"

static PyObject *ErrorObject;

/* ----------------------------------------------------- */

/* Declarations for objects of type SAINVObject */

typedef struct {
  PyObject_HEAD
  /* XXXX Add your own stuff here */
} SAINVObjectobject;

staticforward PyTypeObject SAINVObjecttype;



/* ---------------------------------------------------------------- */

static char SAINVObject_precon__doc__[] = 
""
;

static PyObject *
SAINVObject_precon(self, args)
     SAINVObjectobject *self;
     PyObject *args;
{
  /* parse input parameters */
  if (!PyArg_ParseTuple(args, ""))
    return NULL;


  Py_INCREF(Py_None);
  return Py_None;
}


static struct PyMethodDef SAINVObject_methods[] = {
	{"precon",	(PyCFunction)SAINVObject_precon,	METH_VARARGS,	SAINVObject_precon__doc__},
 
	{NULL,		NULL}		/* sentinel */
};

/* ---------- */


static SAINVObjectobject *
newSAINVObjectobject()
{
	SAINVObjectobject *self;
	
	self = PyObject_NEW(SAINVObjectobject, &SAINVObjecttype);
	if (self == NULL)
		return NULL;
	/* XXXX Add your own initializers here */
	return self;
}


static void
SAINVObject_dealloc(self)
	SAINVObjectobject *self;
{
	/* XXXX Add your own cleanup code here */
	PyMem_DEL(self);
}

static PyObject *
SAINVObject_getattr(self, name)
	SAINVObjectobject *self;
	char *name;
{
	/* XXXX Add your own getattr code here */
	return Py_FindMethod(SAINVObject_methods, (PyObject *)self, name);
}

static char SAINVObjecttype__doc__[] = 
""
;

static PyTypeObject SAINVObjecttype = {
	PyObject_HEAD_INIT(&PyType_Type)
	0,				/*ob_size*/
	"SAINVObject",			/*tp_name*/
	sizeof(SAINVObjectobject),		/*tp_basicsize*/
	0,				/*tp_itemsize*/
	/* methods */
	(destructor)SAINVObject_dealloc,	/*tp_dealloc*/
	(printfunc)0,		/*tp_print*/
	(getattrfunc)SAINVObject_getattr,	/*tp_getattr*/
	(setattrfunc)0,	/*tp_setattr*/
	(cmpfunc)0,		/*tp_compare*/
	(reprfunc)0,		/*tp_repr*/
	0,			/*tp_as_number*/
	0,		/*tp_as_sequence*/
	0,		/*tp_as_mapping*/
	(hashfunc)0,		/*tp_hash*/
	(ternaryfunc)0,		/*tp_call*/
	(reprfunc)0,		/*tp_str*/

	/* Space for future expansion */
	0L,0L,0L,0L,
	SAINVObjecttype__doc__ /* Documentation string */
};

/* End of code for SAINVObject objects */
/* -------------------------------------------------------- */


static char sparslab_sainv__doc__[] =
"sainv(A, droptol=0.1)

parameters:

 A: csr_mat object
 
 droptol: drop tolerance
"
;

static PyObject *
sparslab_sainv(self, args)
     PyObject *self;	/* Not used */
     PyObject *args;
{
  /* parameters */
  PyObject *A;
  double droptol = 0.1;
  /* Fortran parameters for AINVSR2_I */
  int OUTLVL = 1;
  int OUNIT = 6;
  double DRFL = droptol;
  int SIZE_C = 0;
  int DIAG_ONE = -1;		/* output */
  double DIAGTOL = -1.0;
  int DROPTYP = 0;
  double EPS = DBL_EPSILON;	/* constant */
  double MI = -1.0;		/* output */
  int IPAR[600];
  double RPAR[600];
  /* Additional Fortran parameters for AINVSR2 */
  int MSGLVL = 1;
  int MSGUNIT = 6;
  int N;			/* set later */
  int *IA;			/* set later */
  int *JA;			/* set later */
  int *AA;			/* set later */
  IM90;
  JM90;
  AM90;
  AM90;
  int SIZE_P = 0;
  int SIZE_R = 0;
  DIAG9;
  int IMODIF = -1;
  int FILL = 0;			/* output */
  int FILLMAX = 0;		/* output */
  int IFILLMAX = 0;		/* output */
  int GARROW = 0;		/* output */
  int GARCOL = 0;		/* output */
  int INFO = 0;			/* output */
  /* Additional Fortran parameters for AINVSR2_O */
  int PFORMAT = 0;
  int PCLASS = 0;
  
    
  

  if (!PyArg_ParseTuple(args, "O!|d", &CSRMatType, &A, &droptol))
    return NULL;
  
  /* convert A to 1-based indexing */
  
  /* call Fortran routines for SAINV construction */
  F77(AINVSR2_I)(&OUTLVL,&OUNIT,&DRFL,&SIZE_C,&DIAG_ONE,&DIAGTOL,
		 &DROPTYP,&EPS,&MI,IPAR,RPAR);
  F77(AINVSR2)(&MSGLVL,&MSGUNIT,&N,IA,JA,AA,IM90,JM90,AM90,&SIZE_P,
	       &SIZE_C,&SIZE_R,DIAG9,&DIAGTOL,&DRFL,&MI,&DIAG_ONE,
	       &DROPTYP,&IMODIF,&FILL,&FILLMAX,
	       &IFILLMAX,&GARROW,&GARCOL,&INFO)
  F77(AINVSR2_O)(&OUTLVL,&OUNIT,&PFORMAT,&PCLASS,
		 &SIZE_C,&FILLMAX,&IFILLMAX,&IMODIF,&FILL,&GARCOL,&GARROW,&INFO,&IPAR)


  Py_INCREF(Py_None);
  return Py_None;
}

/* List of methods defined in the module */

static struct PyMethodDef sparslab_methods[] = {
  {"sainv", (PyCFunction)sparslab_sainv, METH_VARARGS, sparslab_sainv__doc__},
  {NULL, (PyCFunction)NULL, 0, NULL}		/* sentinel */
};


/* Initialization function for the module (*must* be called initsparslab) */

static char sparslab_module_documentation[] = 
""
;

void
initsparslab()
{
  PyObject *m, *d;
  
  /* Create the module and add the functions */
  m = Py_InitModule4("sparslab", sparslab_methods,
		     sparslab_module_documentation,
		     (PyObject*)NULL,PYTHON_API_VERSION);
  
  /* Add some symbolic constants to the module */
  d = PyModule_GetDict(m);
  ErrorObject = PyString_FromString("sparslab.error");
  PyDict_SetItemString(d, "error", ErrorObject);
  
  /* XXXX Add constants here */

  /* initialize Numeric array module */
  import_array();
  /* initialize spmatrix module */
  import_spmatrix();  
  
  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module sparslab");
}

