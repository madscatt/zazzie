#include "Python.h"
#include "numpy/arrayobject.h"

#include <iostream>

#if CPP_LIB == 1
#include <GV.h>
#include <ScVars.h>
#endif

#if CUDA_LIB == 1
#include <cudaGV.h>
#endif

#if CUDA_DRIVER == 1
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#endif

#define DEBUG 0
#define NGPU 0 ///< @NOTE to ZHL

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static
PyObject *initialize(PyObject *self, PyObject *args)
{
    int Natoms;
    PyObject *py_sascalc_inputs;
    PyObject *py_B_neutron;
    PyObject *py_B_xray;

    ///< pass parameters
	if (!PyArg_ParseTuple(args, "iOOO", &Natoms, &py_sascalc_inputs, &py_B_neutron, &py_B_xray))
    {
        std::cerr<<"Error parsing py_sascalc_inputs from python to cpp extension"<<std::endl;
        return Py_None;
    }

    PyArrayObject *pyArray_B_neutron = (PyArrayObject*)PyArray_ContiguousFromObject(py_B_neutron,PyArray_DOUBLE, 0, 2);
    double *B_neutron = (double*)(pyArray_B_neutron->data);
    if(DEBUG) printf("B_neutron: %f %f %f %f\n",B_neutron[0],B_neutron[1],B_neutron[2],B_neutron[3]);

    PyArrayObject *pyArray_B_xray = (PyArrayObject*)PyArray_ContiguousFromObject(py_B_xray,PyArray_DOUBLE, 0, 3);
    double *B_xray = (double*)(pyArray_B_xray->data);
    if(DEBUG) printf("B_xray: %f %f %f %f\n",B_xray[0],B_xray[1],B_xray[2],B_xray[3]);

    std::string xon(PyBytes_AsString(PyObject_GetAttrString(py_sascalc_inputs, "xon")));
    int neutron_flag=0, xray_flag=0;
    if (xon==std::string("neutron") or xon==std::string("neutron_and_xray")) neutron_flag = 1;
    if (xon==std::string("xray") or xon==std::string("neutron_and_xray")) xray_flag = 1;

    const int Nq = PyLong_AsLong(PyObject_GetAttrString(py_sascalc_inputs, "number_q_values"));
    const double Qmax = PyFloat_AsDouble(PyObject_GetAttrString(py_sascalc_inputs, "q_max"));
    if(DEBUG) printf("Nq: %d, Qmax: %f\n",Nq,Qmax);

    // neutron constrast points data
    std::vector<sascalc::ScVars::neutron_contrast_t> neutron_contrasts_vector;
    if (neutron_flag)
    {
        int Ncontrasts = PyLong_AsLong(PyObject_GetAttrString(py_sascalc_inputs, "number_contrast_points"));
        double D2O_percentage ; //, I0;
        if (Ncontrasts<0) Ncontrasts=0;
        for (int i=0; i<Ncontrasts; ++i)
        {
            D2O_percentage = PyFloat_AsDouble(PyList_GetItem(PyObject_GetAttrString(py_sascalc_inputs, "D2O_percentage_array"), i));
            double I0 = PyFloat_AsDouble(PyList_GetItem(PyObject_GetAttrString(py_sascalc_inputs, "I0_array"), i));
            neutron_contrasts_vector.push_back(std::make_pair(D2O_percentage, I0));
        }
    }

    // exH region data
    sascalc::ScVars::neutron_exH_t neutron_exH;
    std::vector<sascalc::ScVars::neutron_exH_t> neutron_exH_vector;
    const int NexH = PyLong_AsLong(PyObject_GetAttrString(py_sascalc_inputs, "number_exH_regions"));
    std::string exH_basis_string;
    double fraction_exH;
    for (int i=0; i<NexH; ++i)
    {
        //exH_basis_string = PyString_AsString(PyList_GetItem(PyObject_GetAttrString(py_sascalc_inputs, "exH_basis_string_array"),i));
        exH_basis_string = PyBytes_AsString(PyList_GetItem(PyObject_GetAttrString(py_sascalc_inputs, "exH_basis_string_array"),i));
        fraction_exH = PyFloat_AsDouble(PyList_GetItem(PyObject_GetAttrString(py_sascalc_inputs, "fraction_exH_array"),i));
        neutron_exH_vector.push_back(std::make_pair(exH_basis_string, fraction_exH));
    }

    // deuterated region data
    sascalc::ScVars::neutron_deut_t neutron_deut;
    std::vector<sascalc::ScVars::neutron_deut_t> neutron_deut_vector;
    const int Ndeut = PyLong_AsLong(PyObject_GetAttrString(py_sascalc_inputs, "number_deuteration_regions"));
    std::string deut_basis_string;
    double fraction_deut;
    for (int i=0; i<Ndeut; ++i)
    {
        deut_basis_string = PyBytes_AsString(PyList_GetItem(PyObject_GetAttrString(py_sascalc_inputs, "deuterated_basis_string_array"),i));
        fraction_deut = PyFloat_AsDouble(PyList_GetItem(PyObject_GetAttrString(py_sascalc_inputs, "fraction_deuterated_array"),i));
        neutron_deut_vector.push_back(std::make_pair(deut_basis_string, fraction_deut));
    }
    
    // xray constrast points data
    std::vector<sascalc::ScVars::xray_contrast_t> xray_contrasts_vector;
    if (xray_flag)
    {
        int Nxraycontrasts = PyLong_AsLong(PyObject_GetAttrString(py_sascalc_inputs, "xray_number_contrast_points"));
        if (Nxraycontrasts<0) Nxraycontrasts=0;
        for (int i=0; i<Nxraycontrasts; ++i)
        {
            double D2O_percentage = PyFloat_AsDouble(PyList_GetItem(PyObject_GetAttrString(py_sascalc_inputs, "xray_D2O_percentage_array"), i));
            double I0 = PyFloat_AsDouble(PyList_GetItem(PyObject_GetAttrString(py_sascalc_inputs, "xray_I0_array"), i));
            xray_contrasts_vector.push_back(std::make_pair(D2O_percentage, I0));
        }
    }

    // gv data
    const std::string gv_method = PyBytes_AsString(PyObject_GetAttrString(py_sascalc_inputs, "golden_vector_method_option"));
    double gv_parameter;
    if (gv_method == std::string("fixed"))
        gv_parameter = PyLong_AsLong(PyObject_GetAttrString(py_sascalc_inputs, "number_golden_vectors"));
    else if  (gv_method == std::string("converge"))
        gv_parameter = PyFloat_AsDouble(PyObject_GetAttrString(py_sascalc_inputs, "golden_vector_method_converge_tolerance"));

    // construct sascalc::ScVars object
    sascalc::ScVars *p_scvars = new sascalc::ScVars(Nq, Qmax, neutron_contrasts_vector, neutron_exH_vector, neutron_deut_vector, xray_contrasts_vector, gv_method, gv_parameter);
    if(DEBUG) p_scvars->show();

    // construct SasCalc object
    sascalc::SasCalc *p_sascalc;
#if USE_CUDA == 1
    int gpu_device_count;
    cudaError_t cudaErr = cudaGetDeviceCount( &gpu_device_count );
    if (cudaErr == cudaSuccess && gpu_device_count)
    {
        CALL_CUDA(cudaSetDevice(NGPU));
        p_sascalc = new sascalc::cudaGV(Natoms, *p_scvars, B_neutron, B_xray);
        CALL_CUDA(cudaDeviceSynchronize());
        std::cout<<"running Sascalc on GPU"<<std::endl;
    }
    else
    {
        p_sascalc = new sascalc::GV(Natoms, *p_scvars, B_neutron, B_xray);
        std::cout<<"running Sascalc on CPU"<<std::endl;
    }
#elif USE_CPU == 1
    p_sascalc = new sascalc::GV(Natoms, *p_scvars, B_neutron, B_xray);
    std::cout<<"running Sascalc on CPU"<<std::endl;
#endif


    ///< return
    //return PyCObject_FromVoidPtr((void*)p_sascalc, NULL);
    return PyCapsule_New((void*)p_sascalc, NULL, NULL);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// calculator
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static
PyObject *calculate(PyObject *self, PyObject *args)
{
    PyObject * sascalc_object;
    PyObject * py_coordinates;
    int Nframes;

    ///< pass parameters
	if (!PyArg_ParseTuple(args, "OOi", &sascalc_object, &py_coordinates, &Nframes))
    {
        std::cerr<<"Error parsing py_sascalc_inputs from python to cpp extension"<<std::endl;
        Py_INCREF(Py_None);
        return Py_None;
    }

    PyArrayObject * pyArray_coordinates = (PyArrayObject*)PyArray_ContiguousFromObject(py_coordinates,PyArray_DOUBLE, 3,3);
    double *coordinates = (double*)(pyArray_coordinates->data);
    if(DEBUG) printf("Coor: %f %f %f %f\n",coordinates[0],coordinates[1],coordinates[2],coordinates[3]);

    //sascalc::SasCalc *p_sascalc = static_cast<sascalc::GV*>(PyCObject_AsVoidPtr(sascalc_object));
    sascalc::SasCalc *p_sascalc = static_cast<sascalc::GV*>(PyCapsule_GetPointer(sascalc_object, NULL));

    sascalc::ScResults *p_results = p_sascalc -> calculate(coordinates, Nframes);

    ///< return
    //return PyCObject_FromVoidPtr((void*)p_results, NULL);
    return PyCapsule_New((void*)p_results, NULL, NULL);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// save
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static
PyObject *save(PyObject *self, PyObject *args)
{
    PyObject * results_object, *py_name, *py_runname;

    ///< pass parameters
	if (!PyArg_ParseTuple(args, "OOO", &results_object, &py_name, &py_runname))
    {
        std::cerr<<"Error parsing py_sascalc_inputs from python to cpp extension"<<std::endl;
        Py_INCREF(Py_None);
        return Py_None;
    }

    //sascalc::ScResults *p_results = static_cast<sascalc::ScResults*>(PyCObject_AsVoidPtr(results_object));
    sascalc::ScResults *p_results = static_cast<sascalc::ScResults*>(PyCapsule_GetPointer(results_object, NULL));

    const std::string name = PyBytes_AsString(py_name);
    const std::string runname = PyBytes_AsString(py_runname);

    //std::cout<<(*p_results)<<std::endl;
    p_results->save(name, runname);

    ///< return
    //return PyCObject_FromVoidPtr((void*)p_results, NULL);
    return PyCapsule_New((void*)p_results, NULL, NULL);
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
        std::cerr<<"Error parsing py_sascalc_inputs from python to cpp extension"<<std::endl;
        Py_INCREF(Py_None);
        return Py_None;
    }

    //sascalc::SasCalc * p_sascalc = static_cast<sascalc::SasCalc*>(PyCObject_AsVoidPtr(sascalc_object));
    sascalc::SasCalc * p_sascalc = static_cast<sascalc::SasCalc*>(PyCapsule_GetPointer(sascalc_object, NULL));
    delete p_sascalc;

    ///< return
    Py_INCREF(Py_None);
    return Py_None;
}

static PyMethodDef SascalcMethods[] = {
    
    	//{"system",  spam_system, METH_VARARGS,
     	//"Execute a shell command."},
	{ "initialize", initialize, METH_VARARGS },
	{ "calculate", calculate, METH_VARARGS },
	{ "save", save, METH_VARARGS },
	{ "clean", clean, METH_VARARGS },
    	{NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef sascalc_api = {
	PyModuleDef_HEAD_INIT,
    	"sascalc_api",   /* name of module */
    	NULL, /* module documentation, may be NULL */
    	-1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
	SascalcMethods
  
};

PyMODINIT_FUNC
PyInit_sascalc_api(void)
{
	return PyModule_Create(&sascalc_api);
}

/*
PyMODINIT_FUNC
initsascalc_api(){
	PyObject *m;
	m = Py_InitModule("sascalc_api", exampleMethods);
	import_array();
}
*/


