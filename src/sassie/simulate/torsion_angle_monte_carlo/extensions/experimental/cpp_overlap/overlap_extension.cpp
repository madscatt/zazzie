#include "Python.h"
#include "numpy/arrayobject.h"

#include <iostream>
#include <dlfcn.h>

#if CPP_LIB == 1
//#include <overlap.h>
#include "overlap.h"
#endif

#if CUDA_LIB == 1
#include <wrapperCudaKernel.h>
#endif

#if CUDA_DRIVER == 1
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#endif

/*
 * coor is a 2-dimensional numpy array with dimension: natoms * 3
 * r is a 1-dimensional numpy array with dimension: natoms
 * bond_list is a 2-dimensional numpy array with dimension: nbonds * 2
 */

PyObject *overlap(PyObject *self, PyObject *args){

	PyObject *input_coor = NULL ;
	PyObject *input_r = NULL ;
	PyObject *input_bond_list = NULL ;
    int natoms;
    int nbonds;
	double scale;
	
    ///< pass parameters
	if (!PyArg_ParseTuple(args, "OOid|Oi", &input_coor, &input_r, &natoms, &scale, &input_bond_list, &nbonds)) return NULL;

    ///< get coor
    PyArrayObject * py_coor;
    py_coor = (PyArrayObject*)PyArray_ContiguousFromObject(input_coor,PyArray_DOUBLE, 2,2);
    double *coor = (double*)(py_coor->data);
    //printf("test %f %f %f %f\n",coor[0],coor[1],coor[2],coor[3]);

    ///< get radius
    PyArrayObject * py_r;
    py_r = (PyArrayObject*)PyArray_ContiguousFromObject(input_r,PyArray_DOUBLE, 1, 1);
    double *r = (double*)(py_r->data);
    //printf("test %f %f %f %f\n",r[0],r[1],r[2],r[3]);

    ///< get bond_list
    PyArrayObject * py_bond_list=0;
    int *bond_list = 0;
    if (input_bond_list)
    {
        py_bond_list = (PyArrayObject*)PyArray_ContiguousFromObject(input_bond_list,PyArray_INT, 2,2);
        bond_list = (int*)(py_bond_list->data);
        //printf("test %d %d %d %d\n",bond_list[0],bond_list[1],bond_list[2],bond_list[3]);
    }

    // local variables
    int flag_overlap=0 ;
	int i,j;

    ///< build the bond table
    // ----------------->idx_bond_small
    // |------------
    // |*-----------
    // |**----------
    // |***---------
    // |****--------
    // |*****-------
    // |******------
    // |*******-----
    // |********----
    // |*********---
    // |**********--
    // |***********-
    // V
    // idx_bond_big
    // Where '*' is filled and '-' is not filled
    // Note that the diagonal is also not filled
    int size_bond_table=0;
    int * bond_table=0;
    if (py_bond_list)
    {
        size_bond_table = (natoms*(natoms-1))/2;
        bond_table = (int*)malloc(size_bond_table*sizeof(int));
        for (i=0; i<size_bond_table; ++i) bond_table[i] = 0;
        int idx_bond_big,idx_bond_small,tmp;
        for (i=0; i<nbonds; ++i)
        {
            idx_bond_big = bond_list[2*i];
            idx_bond_small = bond_list[2*i+1];
            if (idx_bond_big<idx_bond_small) {tmp=idx_bond_big; idx_bond_big=idx_bond_small; idx_bond_small=tmp;}
            //printf("test bond_big: %d bond_small: %d\n",idx_bond_big,idx_bond_small);
            bond_table[(idx_bond_big*(idx_bond_big-1))/2 + idx_bond_small] = 1;
        }
        /*
        int i_sofar = 0, count=1;
        for (i=0; i<size_bond_table; ++i)
        {
            printf("%d ",bond_table[i]);
            i_sofar ++;
            if (i_sofar%count==0)
            {
                printf("\n");
                count++;
                i_sofar=0;
            }
        }
        */
    }

    ///< some print-outs
//    std::cout<<"natoms: "<<natoms<<std::endl;
//    if (bond_table) std::cout<<"nbonds: "<<nbonds<<std::endl;

#if CUDA_LIB == 1
    ///< Query GPU deivces
    int ngpus = 0;
    cudaError_t cudaErr = cudaGetDeviceCount(&ngpus);
    if (cudaErr != cudaSuccess) ngpus=-1;
    //else std::cout<<"Number of GPUs: "<<ngpus<<std::endl;
    if (ngpus > 0)
    {
        std::cout<<std::endl<<"Checking for overlap on GPU..."<<std::endl;
        double *gpu_coor, *gpu_r;
        int *gpu_bond_table;
        CALL_CUDA(cudaMalloc((void**)&gpu_coor,natoms*3*sizeof(double)));
        CALL_CUDA(cudaMemcpy(gpu_coor,coor,natoms*3*sizeof(double),cudaMemcpyHostToDevice));
        CALL_CUDA(cudaMalloc((void**)&gpu_r,natoms*sizeof(double)));
        CALL_CUDA(cudaMemcpy(gpu_r,r,natoms*sizeof(double),cudaMemcpyHostToDevice));
        CALL_CUDA(cudaMalloc((void**)&gpu_bond_table,size_bond_table*sizeof(int)));
        CALL_CUDA(cudaMemcpy(gpu_bond_table,bond_table,size_bond_table*sizeof(int),cudaMemcpyHostToDevice));
        int *cpu_noverlap, *gpu_noverlap;
        CALL_CUDA(cudaMallocHost((void**)&cpu_noverlap,sizeof(int)));
        cpu_noverlap[0] = 0;
        CALL_CUDA(cudaMalloc((void**)&gpu_noverlap,sizeof(int)));
        CALL_CUDA(cudaMemcpy(gpu_noverlap,cpu_noverlap,sizeof(int),cudaMemcpyHostToDevice));
        wrapper_cudaOverlap(gpu_coor, gpu_r, natoms, scale, gpu_noverlap, gpu_bond_table);
        CALL_CUDA(cudaMemcpy(cpu_noverlap,gpu_noverlap,sizeof(int),cudaMemcpyDeviceToHost));
        if (cpu_noverlap[0]) flag_overlap = 1;
//        if (!flag_overlap) std::cout<<"No overlap found!"<<std::endl;
//        else std::cout<<"Overlap found!"<<std::endl<<"Number of overlaps found: "<<cpu_noverlap[0]<<std::endl;
        CALL_CUDA(cudaFree(gpu_coor));
        CALL_CUDA(cudaFree(gpu_r));
        CALL_CUDA(cudaFree(gpu_bond_table));
        CALL_CUDA(cudaFree(gpu_noverlap));
        CALL_CUDA(cudaFreeHost(cpu_noverlap));
    }
    else
    {
    #ifndef CPP_LIB
        {
            std::cout<<"Only cuda library provided, but no cuda device found!"<<std::endl;
            exit(1);
        }
    #else
        ///< running on CPU
 //       std::cout<<std::endl<<"Checking for overlap on CPU..."<<std::endl;
        flag_overlap = overlap(coor, r, natoms, scale, bond_table);
 //       if (!flag_overlap) std::cout<<"No overlap found!"<<std::endl;
 //       else std::cout<<"Overlap found!"<<std::endl;
    #endif
    }
#else
    ///< running on CPU
 //   std::cout<<std::endl<<"Checking for overlap on CPU..."<<std::endl;
    flag_overlap = overlap(coor, r, natoms, scale, bond_table);
//    if (!flag_overlap) std::cout<<"No overlap found!"<<std::endl;
//    else std::cout<<"Overlap found!"<<std::endl;
#endif


    ///< clean
    free(bond_table);

    ///< return
	return Py_BuildValue("i",flag_overlap) ;
}

static struct PyModuleDef exampleMethods = 
{   
    PyModuleDef_HEAD_INIT,
    "overlap",
    "",
    -1,
    module_methods
};

PyMODINIT_FUNC PyInit_overlap(void)
{
    return PyModule_Create(&overlap) ;
};

//
//  Python 2.X 
//
//static PyMethodDef exampleMethods[] = {
//	{ "overlap", overlap, METH_VARARGS },
//	{ NULL, NULL }
//} ;

//PyMODINIT_FUNC
//initoverlap(){
//	PyObject *m, *m2;
//	m = Py_InitModule("overlap", exampleMethods);
//	import_array();
//}

