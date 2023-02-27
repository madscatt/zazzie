#include <math.h>
#include "Python.h"
#include "numpy/arrayobject.h"

#define IND1(a, i) *((double *)(a->data + i*a->strides[0]))
#define IND2(a, i, j) *((double *)(a->data + i*a->strides[0] + j*a->strides[1]))
#define IND3(a, i, j, k) *((double *)(a->data + i*a->strides[0] + j*a->strides[1] + k*a->strides[2]))

/*
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
*/

PyObject *renorm(PyObject *self, PyObject *args){
	PyObject *input1 ;
	PyArrayObject *work ;
	int i,j,k,nxgp,nygp,nzgp ;
	double val,max ; 

        if (!PyArg_ParseTuple(args, "Oiii", &input1, &nxgp, &nygp, &nzgp)){
		return NULL;
	}	
	work = (PyArrayObject *) PyArray_ContiguousFromObject(input1, PyArray_DOUBLE, 3,3);
	if(work==NULL){
		return NULL; 
	}
	max=0.0 ;
	for(i=0; i<nxgp; i++){
		for(j=0; j<nygp; j++){
			for(k=0; k<nzgp; k++){
				val = IND3(work,i,j,k) ;
				if(val>max)
					max=val ;

			}

		}

	}

	Py_DECREF(work) ;
	return Py_BuildValue("f",max) ;

}

static PyMethodDef renormMethods[] = {
	{ "renorm", renorm, METH_VARARGS },
	{ NULL, NULL }
} ;

void initrenorm(){
	PyObject *m ;
	m = Py_InitModule("renorm", renormMethods);
	import_array();
}

