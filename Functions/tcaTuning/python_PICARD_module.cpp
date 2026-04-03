#include <Python.h>
#include <dace/dace.h>
#include <cmath>
#include <fstream>
#include "numpy/ndarraytypes.h"
#include "numpy/ufuncobject.h"
#include "numpy/npy_3kcompat.h"
#include <PICARD_methods.h>

#define MAXORDER 6
#define MAXVARS 1

using namespace std;

// definition of the python module properties
static PyMethodDef DAgenmethods[] = {
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static PyModuleDef picardmodule = {
    PyModuleDef_HEAD_INIT,
    "picard",
    "Module implementing relevant functions for the implementation of the picard lindelhof expansion method.",
    -1,
    DAgenmethods,
};

// definition of the object methods
static PyMethodDef Picard_methods[] = {
    {"extractcoefficients", (PyCFunction) Picard_extractcoefficients, METH_VARARGS | METH_KEYWORDS,
     "expand in time about desired initial state"
    },
    {NULL}  /* Sentinel */
};

// definition of the pyobject subclass
static PyTypeObject PicardType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "picard.Picard",               /* tp_name */
    sizeof(PicardObject),              /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor) Picard_dealloc,       /* tp_dealloc */
    0,                              /* tp_vectorcall_offset */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_as_async */
    0,                              /* tp_repr */
    0,                              /* tp_as_number */
    0,                              /* tp_as_sequence */
    0,                              /* tp_as_mapping */
    0,                              /* tp_hash */
    0,                              /* tp_call */
    0,                              /* tp_str */
    0,                              /* tp_getattro */
    0,                              /* tp_setattro */
    0,                              /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,             /* tp_flags */
    "Piacrd-Lindelof object",                  /* tp_doc */
    0,                              /* tp_traverse */
    (inquiry) Picard_clear,            /* tp_clear */
    0,                              /* tp_richcompare */
    0,                              /* tp_weaklistoffset */
    0,                              /* tp_iter */
    0,                              /* tp_iternext */
    Picard_methods,                    /* tp_methods */
    0,                              /* tp_members */
    0,                              /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    (initproc) Picard_init,            /* tp_init */
    0,                              /* tp_alloc */
    Picard_new,                        /* tp_new */
};


//////////// Here the methods associated to the object begin 
// definition of method to print to screen polynomial maps:

// destructor
static int
Picard_clear(PicardObject *self)
{
    // delete what is pointed by the address saved in self.state and self.delta_P
    return 0;
}

// free memory
static void
Picard_dealloc(PicardObject *self)
{
    // frees memory of instance as well as the pointer to itself
    //PyObject_GC_UnTrack(self);
    Picard_clear(self);
    Py_TYPE(self)->tp_free((PyObject *) self);
}

// constructor
static int
Picard_init(PicardObject *self, PyObject *args, PyObject *kwds)
{
    return 0;
}

// memory allocator
static PyObject *
Picard_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PicardObject *self;
    self = (PicardObject *) type->tp_alloc(type, 0); // allocate memory pointer for new dagmm type object
    // if memory allocation was successful allocate place for state and covariance polynomial maps
    return (PyObject *) self;
}

void picard_lindelof(DACE::AlgebraicVector<DACE::DA> &x, int ord,  DACE::AlgebraicVector<DACE::DA> (*f)(AlgebraicVector<DACE::DA>,DACE::DA,double), double mu)
//void picard_lindelof(AlgebraicVector<DA> &x, int ord,  AlgebraicVector<DA> (*f)(AlgebraicVector<DA>,DA,Dynamics_Eq*), Dynamics_Eq* context, double scale_fact )
//void picard_lindelof(AlgebraicVector<DA> &x, int ord,  AlgebraicVector<DA> (*f)(AlgebraicVector<DA>,double,Dynamics_Eq*), Dynamics_Eq* context, double scale_fact )
{
		
	DACE::DA t = DACE::DA(1);
	DACE::AlgebraicVector<DACE::DA> x_i(12), x_j(12);
	x_i = x;	
	
	for (int i = 0; i < ord; ++i) {
		// iteration n times
		//x_j = x + f(x_i, t, context).integ(13);
		x_j = x + f(x_i, t, mu).integ(1);
		x_i = x_j;
	}
	
	x=x_i;
}

// function describing Keplerian dynamics
DACE::AlgebraicVector <DACE::DA> KEPLER(DACE::AlgebraicVector<DACE::DA> x ,DACE::DA t, double mu)
{
	DACE::AlgebraicVector <DACE::DA> dx(12);// dummy variable to output first order dynamical equations
    DACE::AlgebraicVector <DACE::DA> r1(3),r2(3); // dummy variables to store position vector and thrust vector;

    // initialization of position and Thrust vector;
    r1[0]=x[0];r1[1]=x[1];r1[2]=x[2];
	r2[0]=x[6];r2[1]=x[7];r2[2]=x[8];
	
    // position vector dynamics (x,y,z)
    dx[0]=x[3];
    dx[1]=x[4];
    dx[2]=x[5];

    // velocity vector dynamics (vx, vy, vz)
    dx[3]=-mu/pow(r1.vnorm(),3.)*x[0];
	dx[4]=-mu/pow(r1.vnorm(),3.)*x[1];
	dx[5]=-mu/pow(r1.vnorm(),3.)*x[2];
	
	dx[6]=x[9];
    dx[7]=x[10];
    dx[8]=x[11];
	
	dx[9]=-mu/pow(r2.vnorm(),3.)*x[6];
	dx[10]=-mu/pow(r2.vnorm(),3.)*x[7];
	dx[11]=-mu/pow(r2.vnorm(),3.)*x[8];
	
	return dx;	
}

// Propagate method
static PyObject *
Picard_extractcoefficients(PicardObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *input=NULL,*input2=NULL;
    PyObject *x10=NULL, *x20=NULL;
    double mu = 398600.4415;

    static char *kwlist[] = {"","","mu","nonlinear_th", "weight_th","N_split", "tolABS", "tolREL", NULL};

    if (!PyArg_ParseTupleAndKeywords(args,kwds, "O!O!|d", (char **)kwlist,&PyArray_Type, &input, &PyArray_Type, &input2, &mu)) return NULL; // parse tuple by order, 1 input is mandatory (generic pyobj), the second optional one can also be specified but it is a numpy array (O!)

    // check first input type, if it is a number it is considered the final propagation time of a thrust free propagatin (int!= float), that is why different checks are implemented
        // if the input is a thrust array  then simply converts it to proper numpy type
    x10 = PyArray_FROM_OTF(input, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY | NPY_F_CONTIGUOUS);
    x20 = PyArray_FROM_OTF(input2, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY | NPY_F_CONTIGUOUS);

    if (x10 == NULL || x20 == NULL) {        //PyArray_DiscardWritebackIfCopy(oarr);
        return NULL;
    }

    double * dataptr_primary = (double*) PyArray_DATA(x10);
    double * dataptr_secondary = (double*) PyArray_DATA(x20);

    DACE::AlgebraicVector<DACE::DA> x_temp(6), x1f(6), x2f(6), d(3),dv(3);

    std::copy(dataptr_primary,dataptr_primary+6,x1f.data());
    std::copy(dataptr_secondary,dataptr_secondary+6,x2f.data());
    x_temp=x1f;
    x_temp<<x2f;

    // decrease reference counters of numpy arrays to avoid memory leaks
    
    picard_lindelof(x_temp, MAXORDER, &KEPLER, mu);

    for (int i=0;i<3;i++) d[i]=x_temp[i]-x_temp[i+6];
    for (int i=3;i<6;i++) dv[i-3]=x_temp[i]-x_temp[i+6];

    DACE::DA r;
	r=d.dot(dv);

    //////////////////////////////////////////////////////////////////////
    std::vector<unsigned int> ord0(1),ord1(1),ord2(1),ord3(1),ord4(1);
    ord0[0]=0; ord1[0]=1; ord2[0]=2; ord3[0]=3; ord4[0]=4;

    Py_DECREF(x10);
    Py_DECREF(x20);

    return Py_BuildValue("ddddd",r.getCoefficient(ord0),r.getCoefficient(ord1),r.getCoefficient(ord2),r.getCoefficient(ord3),r.getCoefficient(ord4));
}

// module initializer function
PyMODINIT_FUNC PyInit_PICARD(void)
{
    // module initializer needs to initialize DA to custom order and number of variable 
    DA::init(MAXORDER, MAXVARS);
    DA::setEps(1e-40);
    // DA::setTO(2);

    // important to have functioning numpy arrays in c++
    import_array();

    // declaration of pyobject instance
    PyObject *m;
    if (PyType_Ready(&PicardType) < 0)
        return NULL;
    

    m = PyModule_Create(&picardmodule);
    Py_INCREF(&PicardType);
    PyModule_AddObject(m, "picardlindelof", (PyObject *)&PicardType);

    return m;

}