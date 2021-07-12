bool flag_Debye = 0;
//bool flag_Debye = 1;

#include "Python.h"
#include "numpy/arrayobject.h"

#include <iostream>
#include <dlfcn.h>

#if CPP_LIB == 1
#include <GV.h>
#include <Debye.h>
#endif

#if CUDA_LIB == 1
#include <cudaGV.h>
#include <cudaDebye.h>
#endif

#if CUDA_DRIVER == 1
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#endif

#define NGPU 0 ///< @NOTE to ZHL

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static
PyObject *initialize(PyObject *self, PyObject *args)
{
    PyObject * sascalc_inputs;

    ///< pass parameters
	if (!PyArg_ParseTuple(args, "O", &sascalc_inputs))
    {
        std::cerr<<"Error parsing sascalc_inputs from python to cpp extension"<<std::endl;
        return Py_None;
    }

    PyObject* py_coordinates = PyObject_GetAttrString(sascalc_inputs, "coor");
    PyArrayObject * pyArray_coordinates;
    pyArray_coordinates = (PyArrayObject*)PyArray_ContiguousFromObject(py_coordinates,PyArray_DOUBLE, 3,3);
    double *coordinates = (double*)(pyArray_coordinates->data);
    //printf("Coor: %f %f %f %f\n",coordinates[0],coordinates[1],coordinates[2],coordinates[3]);

    PyObject* py_Q = PyObject_GetAttrString(sascalc_inputs, "Q");
    PyArrayObject * pyArray_Q;
    pyArray_Q = (PyArrayObject*)PyArray_ContiguousFromObject(py_Q,PyArray_DOUBLE, 1,1);
    double *Q = (double*)(pyArray_Q->data);
    //printf("Q: %f %f %f %f\n",Q[0],Q[1],Q[2],Q[3]);
    
    PyObject* py_R = PyObject_GetAttrString(sascalc_inputs, "R");
    PyArrayObject * pyArray_R;
    pyArray_R = (PyArrayObject*)PyArray_ContiguousFromObject(py_R,PyArray_DOUBLE, 1,1);
    double *R = (double*)(pyArray_R->data);
    //printf("R: %f %f %f %f\n",R[0],R[1],R[2],R[3]);

    //int Nframes = PyArray_DIMS(pyArray_coordinates)[0];
    int Nframes = PyInt_AsLong(PyObject_GetAttrString(sascalc_inputs, "frames_per_batch"));
    int Natoms = PyArray_DIMS(pyArray_coordinates)[2];
    int Nq = PyArray_DIMS(pyArray_Q)[0];
    int Nr = PyArray_DIMS(pyArray_R)[0];
    double dr = R[1]; //< @NOTE to developer: assuming R is a mesh grid starting from 0

    sascalc::SasCalc * p_gv;
#if USE_CUDA == 1
    int gpu_device_count;
    cudaError_t cudaErr = cudaGetDeviceCount( &gpu_device_count );
    if (cudaErr == cudaSuccess && gpu_device_count)
    {
        CALL_CUDA(cudaSetDevice(NGPU));
        if (flag_Debye)
        {
            p_gv = new sascalc::cudaDebye(Natoms, Nframes, Nq, Nr, dr, coordinates, Q);
            std::cout<<"running Debye using CUDA"<<std::endl;
        }
        else
        {
            p_gv = new sascalc::cudaGV(Natoms, Nframes, Nq, Nr, dr, coordinates, Q);
            std::cout<<"running GV using CUDA"<<std::endl;
        }
        CALL_CUDA(cudaDeviceSynchronize());
    }
    else
    {
        if (flag_Debye)
        {
            p_gv = new sascalc::Debye(Natoms, Nframes, Nq, Nr, dr, coordinates, Q);
            std::cout<<"running Debye on CPU"<<std::endl;
        }
        else
        {
            p_gv = new sascalc::GV(Natoms, Nframes, Nq, Nr, dr, coordinates, Q);
            std::cout<<"running GV on CPU"<<std::endl;
        }
    }
#elif USE_CPU == 1
    if (flag_Debye)
    {
        p_gv = new sascalc::Debye(Natoms, Nframes, Nq, Nr, dr, coordinates, Q);
        std::cout<<"running Debye on CPU"<<std::endl;
    }
    else
    {
        p_gv = new sascalc::GV(Natoms, Nframes, Nq, Nr, dr, coordinates, Q);
        std::cout<<"running GV on CPU"<<std::endl;
    }
#endif

    ///< return
    return PyCObject_FromVoidPtr((void*)p_gv, NULL);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static
PyObject *
batch_load(PyObject *self, PyObject *args)
{
    PyObject * sascalc_object;
    PyObject * py_start;
    PyObject * py_extend;

    ///< pass parameters
	if (!PyArg_ParseTuple(args, "OOO", &sascalc_object, &py_start, &py_extend))
    {
        std::cerr<<"Error parsing sascalc_inputs from python to cpp extension"<<std::endl;
        Py_INCREF(Py_None);
        return Py_None;
    }

    const int start = PyInt_AsLong(py_start);
    const int extend = PyInt_AsLong(py_extend);

    sascalc::SasCalc * p_gv = static_cast<sascalc::SasCalc*>(PyCObject_AsVoidPtr(sascalc_object));
    p_gv->batch_load(start, extend);

    ///< return
    Py_INCREF(Py_None);
    return Py_None;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// a help function for calculate
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static
int tolerance_to_ngv(double tolerance)
{
    if (tolerance>=0.05) return 11;
    if (tolerance<0.05 && tolerance>=0.04) return 13;
    if (tolerance<0.04 && tolerance>=0.03) return 15;
    if (tolerance<0.03 && tolerance>=0.02) return 19;
    if (tolerance<0.02 && tolerance>=0.01) return 31;
    if (tolerance<0.01) return 35;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static
PyObject *calculate(PyObject *self, PyObject *args)
{
    PyObject * sascalc_object;
    PyObject * sascalc_inputs;
    PyObject * sascalc_results;
    PyObject * py_frame;

    ///< pass parameters
	if (!PyArg_ParseTuple(args, "OOOO", &sascalc_object, &py_frame, &sascalc_inputs, &sascalc_results))
    {
        std::cerr<<"Error parsing sascalc_inputs from python to cpp extension"<<std::endl;
        Py_INCREF(Py_None);
        return Py_None;
    }

    PyObject* py_Pr = PyObject_GetAttrString(sascalc_results, "Pr");
    PyArrayObject * pyArray_Pr;
    pyArray_Pr = (PyArrayObject*)PyArray_ContiguousFromObject(py_Pr,PyArray_DOUBLE, 1,1);
    double *Pr_out = (double*)(pyArray_Pr->data);
    int Nr = PyArray_DIMS(pyArray_Pr)[0];
    
    std::string xon(PyString_AsString(PyObject_GetAttrString(sascalc_inputs, "xon")));
    int neutron_flag=0, xray_flag=0;
    if (xon==std::string("neutron") or xon==std::string("neutron_and_xray")) neutron_flag = 1;
    if (xon==std::string("xray") or xon==std::string("neutron_and_xray")) xray_flag = 1;

    const double *B_neutron_array;
    PyObject* py_B_neutron_array = PyObject_GetAttrString(sascalc_inputs, "B_neutron_array");
    PyArrayObject * pyArray_B_neutron_array;
    pyArray_B_neutron_array = (PyArrayObject*)PyArray_ContiguousFromObject(py_B_neutron_array,PyArray_DOUBLE, 3,3);
    B_neutron_array = (double*)(pyArray_B_neutron_array->data);
    const int Ncontrast_neutron  = PyArray_DIMS(pyArray_B_neutron_array)[1];
    const int Natoms = PyArray_DIMS(pyArray_B_neutron_array)[2];

    const double *B_xray_array;
    PyObject* py_B_xray_array = PyObject_GetAttrString(sascalc_inputs, "B_xray_array");
    PyArrayObject * pyArray_B_xray_array;
    pyArray_B_xray_array = (PyArrayObject*)PyArray_ContiguousFromObject(py_B_xray_array,PyArray_DOUBLE, 4,4);
    B_xray_array = (double*)(pyArray_B_xray_array->data);
    const int Ncontrast_xray  = PyArray_DIMS(pyArray_B_xray_array)[1];

    double *Iq_neutron_array;
    PyObject* py_Iq_neutron_array = PyObject_GetAttrString(sascalc_results, "Iq_neutron_array");
    PyArrayObject * pyArray_Iq_neutron_array;
    pyArray_Iq_neutron_array = (PyArrayObject*)PyArray_ContiguousFromObject(py_Iq_neutron_array,PyArray_DOUBLE, 3,3);
    Iq_neutron_array = (double*)(pyArray_Iq_neutron_array->data);
    const int Nq = PyArray_DIMS(pyArray_Iq_neutron_array)[2];

    double *Iq_xray_array;
    PyObject* py_Iq_xray_array = PyObject_GetAttrString(sascalc_results, "Iq_xray_array");
    PyArrayObject * pyArray_Iq_xray_array;
    pyArray_Iq_xray_array = (PyArrayObject*)PyArray_ContiguousFromObject(py_Iq_xray_array,PyArray_DOUBLE, 3,3);
    Iq_xray_array = (double*)(pyArray_Iq_xray_array->data);

    int32_t *Ngv_converged_neutron_array;
    PyObject* py_Ngv_converged_neutron_array = PyObject_GetAttrString(sascalc_results, "Ngv_converged_neutron_array");
    PyArrayObject * pyArray_Ngv_converged_neutron_array;
    pyArray_Ngv_converged_neutron_array = (PyArrayObject*)PyArray_ContiguousFromObject(py_Ngv_converged_neutron_array,PyArray_INT32, 2,2);
    Ngv_converged_neutron_array = (int32_t*)(pyArray_Ngv_converged_neutron_array->data);

    int32_t *Ngv_converged_xray_array;
    PyObject* py_Ngv_converged_xray_array = PyObject_GetAttrString(sascalc_results, "Ngv_converged_xray_array");
    PyArrayObject * pyArray_Ngv_converged_xray_array;
    pyArray_Ngv_converged_xray_array = (PyArrayObject*)PyArray_ContiguousFromObject(py_Ngv_converged_xray_array,PyArray_INT32, 2,2);
    Ngv_converged_xray_array = (int32_t*)(pyArray_Ngv_converged_xray_array->data);

    const int frame = PyInt_AsLong(py_frame);

    sascalc::SasCalc * p_gv = static_cast<sascalc::SasCalc*>(PyCObject_AsVoidPtr(sascalc_object));

    int Ngv;
    double * Iq, *Iq_real, *Iq_imag;
    double * Pr;
    int i, idx, idx_contrast, offset, offset_real, offset_imag;
    const double *B;
    if (Ncontrast_xray>0) xray_flag = 1;
    else xray_flag = 0;
    std::string option(PyString_AsString(PyObject_GetAttrString(sascalc_inputs, "golden_vector_method_option")));
    if (option==std::string("fixed"))
    {
        Ngv = PyInt_AsLong(PyObject_GetAttrString(sascalc_inputs, "number_golden_vectors"));
        // neutron
        for (idx_contrast=0; idx_contrast!=Ncontrast_neutron; ++idx_contrast)
        {
            for (idx=0; idx!=3; ++idx)
            {
                B = &B_neutron_array[idx*(Ncontrast_neutron*Natoms) + idx_contrast*Natoms];
                if (!flag_Debye) static_cast<sascalc::GV*>(p_gv)->updateGV(Ngv);
                p_gv->calculate(B, frame);
                Iq = p_gv->Iq();
                if (!flag_Debye)
                {
                    Iq_real = static_cast<sascalc::GV*>(p_gv)->Iq_real();
                    Iq_imag = static_cast<sascalc::GV*>(p_gv)->Iq_imag();
                }
                else
                {
                    Iq_real = (double *)calloc(Nq, sizeof(double));
                    Iq_imag = (double *)calloc(Nq, sizeof(double));
                }
                offset = idx*(Ncontrast_neutron*Nq) + idx_contrast*Nq;
                offset_real = (idx+3)*(Ncontrast_neutron*Nq) + idx_contrast*Nq;
                offset_imag = (idx+6)*(Ncontrast_neutron*Nq) + idx_contrast*Nq;
                for (i=0; i<Nq; ++i)
                {
                    Iq_neutron_array[offset+i]=Iq[i];
                    Iq_neutron_array[offset_real+i]=Iq_real[i];
                    Iq_neutron_array[offset_imag+i]=Iq_imag[i];
                }
            }
        }
        // xray
        for (idx_contrast=0; idx_contrast!=Ncontrast_xray; ++idx_contrast)
        {
            for (idx=0; idx!=3; ++idx)
            {
                B = &B_xray_array[idx*(Ncontrast_xray*Nq*Natoms) + idx_contrast*Nq*Natoms];
                if (!flag_Debye) static_cast<sascalc::GV*>(p_gv)->updateGV(Ngv);
                p_gv->calculate(B, frame, xray_flag);
                Iq = p_gv->Iq();
                if (!flag_Debye)
                {
                    Iq_real = static_cast<sascalc::GV*>(p_gv)->Iq_real();
                    Iq_imag = static_cast<sascalc::GV*>(p_gv)->Iq_imag();
                }
                else
                {
                    Iq_real = (double *)calloc(Nq, sizeof(double));
                    Iq_imag = (double *)calloc(Nq, sizeof(double));
                }
                offset = idx*(Ncontrast_xray*Nq) + idx_contrast*Nq;
                offset_real = (idx+3)*(Ncontrast_xray*Nq) + idx_contrast*Nq;
                offset_imag = (idx+6)*(Ncontrast_xray*Nq) + idx_contrast*Nq;
                for (i=0; i<Nq; ++i)
                {
                    Iq_xray_array[offset+i]=Iq[i];
                    Iq_xray_array[offset_real+i]=Iq_real[i];
                    Iq_xray_array[offset_imag+i]=Iq_imag[i];
                }
            }
        }
        //p_gv->calculate_pr(frame);
        Pr = p_gv->Pr();
        for (i=0; i<Nr; ++i) Pr_out[i]=Pr[i];
    }
    else if (option==std::string("converge"))
    {
        bool flag_converged=false;
        double * Iq_run_ave = (double*)malloc(Nq*sizeof(double));
        double * Iq_run_ave_previous = (double*)malloc(Nq*sizeof(double));
        double err;
        const double tolerance_gv = PyFloat_AsDouble(PyObject_GetAttrString(sascalc_inputs, "golden_vector_method_converge_tolerance")); 
        const int Ngv_start = tolerance_to_ngv(tolerance_gv);
        double * Iq_result;
        // neutron
        for (idx_contrast=0; idx_contrast!=Ncontrast_neutron; ++idx_contrast)
        {
            for (idx=0; idx!=3; ++idx)
            {
                int count=1;
                flag_converged = false;
                for (Ngv=Ngv_start; !flag_converged && count<=100; Ngv+=2) ///< @NOTE to developers: restrict the number of trials to be 100
                {
                    B = &B_neutron_array[idx*(Ncontrast_neutron*Natoms) + idx_contrast*Natoms];
                    if (!flag_Debye) static_cast<sascalc::GV*>(p_gv)->updateGV(Ngv);
                    p_gv->calculate(B, frame);
                    Iq_result = p_gv->Iq();
                    if (!flag_Debye)
                    {
                        Iq_real = static_cast<sascalc::GV*>(p_gv)->Iq_real();
                        Iq_imag = static_cast<sascalc::GV*>(p_gv)->Iq_imag();
                    }
                    else
                    {
                        Iq_real = (double *)calloc(Nq, sizeof(double));
                        Iq_imag = (double *)calloc(Nq, sizeof(double));
                    }
                    if (count==1) for (i=0; i<Nq; ++i) Iq_run_ave[i]=Iq_result[i];
                    else
                    {
                        err = 0.0;
                        for (i=0; i<Nq; ++i)
                        {
                            Iq_run_ave_previous[i]=Iq_run_ave[i];
                            Iq_run_ave[i] = (Iq_run_ave[i]*(count-1)+Iq_result[i])/count;
                            err += fabs(Iq_run_ave[i]-Iq_run_ave_previous[i])/Iq_run_ave_previous[i];
                        }
                        err /= Nq;
                        if (err<=tolerance_gv) flag_converged = true;
                    }
                    ++count;
                }
                Ngv -= 2;
                Ngv_converged_neutron_array[idx*Ncontrast_neutron+idx_contrast] = Ngv;
                offset = idx*(Ncontrast_neutron*Nq) + idx_contrast*Nq;
                offset_real = (idx+3)*(Ncontrast_neutron*Nq) + idx_contrast*Nq;
                offset_imag = (idx+6)*(Ncontrast_neutron*Nq) + idx_contrast*Nq;
                for (i=0; i<Nq; ++i)
                {
                    Iq_neutron_array[offset+i]=Iq_result[i];
                    Iq_neutron_array[offset_real+i]=Iq_real[i];
                    Iq_neutron_array[offset_imag+i]=Iq_imag[i];
                }
            }
        }
        // xray
        for (idx_contrast=0; idx_contrast!=Ncontrast_xray; ++idx_contrast)
        {
            for (idx=0; idx!=3; ++idx)
            {
                int count=1;
                flag_converged = false;
                for (Ngv=Ngv_start; !flag_converged && count<=100; Ngv+=2) ///< @NOTE to developers: restrict the number of trials to be 100
                {
                    B = &B_xray_array[idx*(Ncontrast_xray*Nq*Natoms) + idx_contrast*Nq*Natoms];
                    if (!flag_Debye) static_cast<sascalc::GV*>(p_gv)->updateGV(Ngv);
                    p_gv->calculate(B, frame, xray_flag);
                    Iq_result = p_gv->Iq();
                    if (!flag_Debye)
                    {
                        Iq_real = static_cast<sascalc::GV*>(p_gv)->Iq_real();
                        Iq_imag = static_cast<sascalc::GV*>(p_gv)->Iq_imag();
                    }
                    else
                    {
                        Iq_real = (double *)calloc(Nq, sizeof(double));
                        Iq_imag = (double *)calloc(Nq, sizeof(double));
                    }
                    if (count==1) for (i=0; i<Nq; ++i) Iq_run_ave[i]=Iq_result[i];
                    else
                    {
                        err = 0.0;
                        for (i=0; i<Nq; ++i)
                        {
                            Iq_run_ave_previous[i]=Iq_run_ave[i];
                            Iq_run_ave[i] = (Iq_run_ave[i]*(count-1)+Iq_result[i])/count;
                            err += fabs(Iq_run_ave[i]-Iq_run_ave_previous[i])/Iq_run_ave_previous[i];
                        }
                        err /= Nq;
                        if (err<=tolerance_gv) flag_converged = true;
                    }
                    ++count;
                }
                Ngv -= 2;
                Ngv_converged_xray_array[idx*Ncontrast_xray+idx_contrast] = Ngv;
                offset = idx*(Ncontrast_xray*Nq) + idx_contrast*Nq;
                offset_real = (idx+3)*(Ncontrast_xray*Nq) + idx_contrast*Nq;
                offset_imag = (idx+6)*(Ncontrast_xray*Nq) + idx_contrast*Nq;
                for (i=0; i<Nq; ++i)
                {
                    Iq_xray_array[offset+i]=Iq_result[i];
                    Iq_xray_array[offset_real+i]=Iq_real[i];
                    Iq_xray_array[offset_imag+i]=Iq_imag[i];
                }
            }
        }
        //p_gv->calculate_pr(frame);
        Pr = p_gv->Pr();
        for (i=0; i<Nr; ++i) Pr_out[i]=Pr[i];
        PyObject_SetAttrString(sascalc_results, "starting_number_golden_vectors", PyInt_FromLong(Ngv_start));
        PyObject_SetAttrString(sascalc_results, "converged_number_golden_vectors", PyInt_FromLong(Ngv));
        free(Iq_run_ave);
        free(Iq_run_ave_previous);
    }

    ///< return
    Py_INCREF(Py_None);
    return Py_None;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static
PyObject *clean(PyObject *self, PyObject *args)
{
    PyObject * sascalc_object;

    ///< pass parameters
	if (!PyArg_ParseTuple(args, "O", &sascalc_object))
    {
        std::cerr<<"Error parsing sascalc_inputs from python to cpp extension"<<std::endl;
        Py_INCREF(Py_None);
        return Py_None;
    }

    sascalc::SasCalc * p_gv = static_cast<sascalc::SasCalc*>(PyCObject_AsVoidPtr(sascalc_object));
    delete p_gv;

    ///< return
    Py_INCREF(Py_None);
    return Py_None;
}



static PyMethodDef exampleMethods[] = {
	{ "initialize", initialize, METH_VARARGS },
	{ "batch_load", batch_load, METH_VARARGS },
	{ "calculate", calculate, METH_VARARGS },
	{ "clean", clean, METH_VARARGS },
	{ NULL, NULL }
} ;

PyMODINIT_FUNC
initsascalc_api(){
	PyObject *m;
	m = Py_InitModule("sascalc_api", exampleMethods);
	import_array();
}

