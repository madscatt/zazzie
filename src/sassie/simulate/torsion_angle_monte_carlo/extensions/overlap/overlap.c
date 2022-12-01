#include <math.h>
#include "Python.h"
#include "numpy/arrayobject.h"
//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

/*
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY;
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
*/

PyObject *overlap(PyObject *self, PyObject *args){
	PyObject *cut ;
	PyArrayObject *array = NULL ;
	double x1,y1,z1,x2,y2,z2,sdist,dist ;
	float lcut ;
	int i,j,natoms,check;

	if (!PyArg_ParseTuple(args, "O!O", &PyArray_Type, &array, &cut))
		return NULL;
	lcut=PyFloat_AS_DOUBLE(cut) ;
	natoms = array->dimensions[0];
	check=0 ;
	for (i=0; i< natoms-1 ; i++){
		x1=*(float *)(array->data + i*array->strides[0]+0*array->strides[1]) ;
		y1=*(float *)(array->data + i*array->strides[0]+(0+1)*array->strides[1]) ;
		z1=*(float *)(array->data + i*array->strides[0]+(0+2)*array->strides[1]) ;
		for (j=i+1; j<natoms ; j++){
			x2=*(float *)(array->data + j*array->strides[0]+0*array->strides[1]) ;
			y2=*(float *)(array->data + j*array->strides[0]+(0+1)*array->strides[1]) ;
			z2=*(float *)(array->data + j*array->strides[0]+(0+2)*array->strides[1]) ;
			sdist=((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1)) ;
			dist=sqrt(sdist) ;
			if(dist<lcut){
				check=1 ;
				//printf("SCH %i\t%i\t%f\n",i,j,dist);

				break ;
			}
		}
		if(check==1) {
			break ;
		}

	}

	return Py_BuildValue("i",check) ;
}

PyObject *moloverlap(PyObject *self, PyObject *args){
	PyObject *cut ;
	PyArrayObject *array_1 = NULL, *array_2 = NULL, *interres=NULL;
	double x1,y1,z1,x2,y2,z2,sdist,dist ;
	float lcut ;
	int i,j,k,natoms_1,natoms_2,resid_1,resid_2,nex,qex,check;

	if (!PyArg_ParseTuple(args, "O!O!O|O!", &PyArray_Type, &array_1,&PyArray_Type, &array_2, &cut, &PyArray_Type, &interres))
		return NULL;
	lcut=PyFloat_AS_DOUBLE(cut) ;
	natoms_1 = array_1->dimensions[0];
	natoms_2 = array_2->dimensions[0];
    nex = interres->dimensions[0];
    //printf("ZHL %3d\n",nex);
	check=0 ;

	for (i=0; i< natoms_1; i++){
		x1=*(float *)(array_1->data + i*array_1->strides[0]+0*array_1->strides[1]) ;
		y1=*(float *)(array_1->data + i*array_1->strides[0]+(0+1)*array_1->strides[1]) ;
		z1=*(float *)(array_1->data + i*array_1->strides[0]+(0+2)*array_1->strides[1]) ;
		for (j=0; j<natoms_2; j++){
			x2=*(float *)(array_2->data + j*array_2->strides[0]+0*array_2->strides[1]) ;
			y2=*(float *)(array_2->data + j*array_2->strides[0]+(0+1)*array_2->strides[1]) ;
			z2=*(float *)(array_2->data + j*array_2->strides[0]+(0+2)*array_2->strides[1]) ;
			qex=0;
			for (k=0; k<nex; k++)
            {
                resid_1 = *(int*)(interres->data + k*interres->strides[0] + 0*interres->strides[1]);
                resid_2 = *(int*)(interres->data + k*interres->strides[0] + 1*interres->strides[1]);
				if (i==resid_1 && j==resid_2)
                {
					qex=1;
                    break;
                }
            }
            if (qex)
            {
                continue;
            }
            sdist=((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1)) ;
			dist=sqrt(sdist) ;
            //printf("ZHL %8.3f\n",dist);
			if(dist<lcut){
				check=1 ;
				break ;
			}
		}
		if(check==1) {
			break ;
		}

	}

	return Py_BuildValue("i",check) ;
}

static PyMethodDef module_methods[] = {
	{ "overlap", overlap, METH_VARARGS },
	{ "moloverlap", moloverlap, METH_VARARGS },
	{ NULL, NULL }
} ;

/* new for python 3.X */

static struct PyModuleDef c_overlap =
{
    PyModuleDef_HEAD_INIT,
    "c_overlap",
    "",
    -1,
    module_methods
};

PyMODINIT_FUNC PyInit_ooverlap(void)
{
    return PyModule_Create(&c_overlap) ;
}

/* worked for python 2.X */

/*
 
void initooverlap(){
	PyObject *m, *m2;
	m = Py_InitModule("ooverlap", exampleMethods);
	m2 = Py_InitModule("moloverlap", exampleMethods);
	import_array();
}
*/


