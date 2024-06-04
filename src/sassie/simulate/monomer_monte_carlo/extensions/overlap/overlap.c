#include <math.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "Python.h"
#include "numpy/arrayobject.h"

PyObject *overlap(PyObject *self, PyObject *args){
    PyObject *cut ;
    PyArrayObject *array = NULL ;
    double x1,y1,z1,x2,y2,z2,sdist,dist ;
    float lcut ;
    int i,j,natoms,check;
    
    if (!PyArg_ParseTuple(args, "O!O", &PyArray_Type, &array, &cut))
        return NULL;
    if (!PyFloat_Check(cut)) {
        PyErr_SetString(PyExc_TypeError, "cut must be a float");
        return NULL;
    }
    lcut=PyFloat_AS_DOUBLE(cut) ;
    if (PyArray_NDIM(array) != 2 || PyArray_DIM(array, 1) != 3) {
        PyErr_SetString(PyExc_TypeError, "array must be 2D with 3 columns");
        return NULL;
    }
    if (PyArray_TYPE(array) != NPY_FLOAT) {
        PyErr_SetString(PyExc_TypeError, "array must contain floats");
        return NULL;
    }
    natoms = PyArray_DIM(array, 0);
    check=0 ;
    for (i=0; i< natoms-1 ; i++){
        x1=*(float *)(PyArray_GETPTR2(array, i, 0)) ;
        y1=*(float *)(PyArray_GETPTR2(array, i, 1)) ;
        z1=*(float *)(PyArray_GETPTR2(array, i, 2)) ;
        for (j=i+1; j<natoms ; j++){
            x2=*(float *)(PyArray_GETPTR2(array, j, 0)) ;
            y2=*(float *)(PyArray_GETPTR2(array, j, 1)) ;
            z2=*(float *)(PyArray_GETPTR2(array, j, 2)) ;
            sdist=((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1)) ;
            dist=sqrt(sdist) ;


            if(dist<lcut){
				check=1 ;
				printf(" i = %d\tj = %d\tdist = %f\n", i,j,dist) ;
				printf(" x1 = %f\t y1 = %f\t z1 = %f\n", x1,y1,z1) ;
				printf(" x2 = %f\t y2 = %f\t z2 = %f\n", x2,y2,z2) ;
				break ;
			}
		}
		if(check==1) {
			break ;
		}

	}

	return Py_BuildValue("i",check) ;
}



static PyMethodDef exampleMethods[] = {
    { "overlap", overlap, METH_VARARGS, "Calculate overlap" },
    { NULL, NULL, 0, NULL }
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "overlap",
    NULL,
    -1,
    exampleMethods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit_overlap(void) {
    import_array();
    return PyModule_Create(&moduledef);
}
