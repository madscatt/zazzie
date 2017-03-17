#include <math.h>
#include "Python.h"
#include "numpy/arrayobject.h"

#define IND1(a, i) *((float *)(a->data + i*a->strides[0]))
#define IND2(a, i, j) *((float *)(a->data + i*a->strides[0] + j*a->strides[1]))
#define IND3(a, i, j, k) *((float *)(a->data + i*a->strides[0] + j*a->strides[1] + k*a->strides[2]))

/*

    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
*/

PyObject *cube(PyObject *self, PyObject *args){
	PyObject *input1, *input2, *input3, *input4 ;
	PyArrayObject *carb, *work ;
	float *lcarb,***gbin ;
	int i,j,k,c,ix,looking,nxgp,nygp,nzgp,tnc,ncarb ;
	float gs,xc,yc,zc,lxx,lxm,lyx,lym,lzx,lzm ;
	float xmin,ymin,zmin,weight ;

/*	printf("am I even here?\n\n");
	fflush(stdout) ;
*/
        if (!PyArg_ParseTuple(args, "OOfffffiii", &input1, &input2, &weight, &xmin, &ymin, &zmin, &gs, &nxgp, &nygp, &nzgp)){
		printf("parse tuple is null\n\n") ;
		fflush(stdout) ;
		
		return NULL;
	}
	carb = (PyArrayObject *) PyArray_ContiguousFromObject(input1, PyArray_FLOAT, 2,2);
	if(carb==NULL){
		printf("carb is null\n\n") ;
		fflush(stdin) ;
		return NULL;
	}
	work = (PyArrayObject *) PyArray_ContiguousFromObject(input2, PyArray_FLOAT, 3,3);
	if(work==NULL){
		printf("work is null\n\n") ;
		fflush(stdin) ;
		return NULL;
	}
 	ncarb = carb->dimensions[0] ;
	tnc=0 ;
	for(c=0; c<ncarb; c++){
		xc=IND2(carb,c,0) ;
		yc=IND2(carb,c,1) ;
		zc=IND2(carb,c,2) ;
  		looking=1 ;
		while((looking==1)){
			for(i=0; i<nxgp; i++){
				lxx=xmin+((i+1)*gs) ;
				lxm=xmin+(i*gs) ;
				if(xc<lxx && xc>=lxm){
					for(j=0; j<nygp; j++){
						lyx=ymin+((j+1)*gs) ;
						lym=ymin+(j*gs) ;
						if(yc<lyx && yc>=lym){
							for(k=0; k<nzgp; k++){
								lzx=zmin+((k+1)*gs) ;
								lzm=zmin+(k*gs) ;
								if(zc<lzx && zc>=lzm){
									IND3(work,i,j,k) = IND3(work,i,j,k)+(1.0*weight) ;
									tnc++ ;
									looking=0 ;	
								}
							}
						}
					}
				}
			

			}
			if(looking==1){
				printf("carbon = %i\n",c) ;
				printf("NO GRIDPOINT FOUND!!\n") ;
				printf("%lf\t%lf\t%lf\n",xc,yc,zc) ;
				return NULL ;
				}

		}

	}

	Py_DECREF(carb) ;

	return PyArray_Return(work) ;
}

static PyMethodDef cubeMethods[] = {
	{ "cube", cube, METH_VARARGS },
	{ NULL, NULL }
} ;

void initcube(){
	PyObject *m ;
	m = Py_InitModule("cube", cubeMethods);
	import_array();
}

