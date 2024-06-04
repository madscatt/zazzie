#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>

static PyObject* pairs(PyObject* self, PyObject* args) {
    PyObject *avdw_list;
    PyArrayObject *cutoffarray_array;
    double r1, r2;
    npy_intp natoms, npairs;

    if (!PyArg_ParseTuple(args, "O!O!", &PyList_Type, &avdw_list, &PyArray_Type, &cutoffarray_array)) {
        return NULL;
    }

    natoms = PyList_Size(avdw_list);
    npairs = PyArray_DIM(cutoffarray_array, 0);

    for (npy_intp i = 0; i < natoms - 1; i++) {
        r1 = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(avdw_list, i), 1));
        for (npy_intp j = i + 1; j < natoms; j++) {
            r2 = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(avdw_list, j), 1));
            *(double *)PyArray_GETPTR1(cutoffarray_array, i + j) = r1 + r2;
        }
    }

    Py_RETURN_NONE;
}

static PyMethodDef module_methods[] = {
    {"pairs", pairs, METH_VARARGS, "Calculate pairs"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "pairs",
    NULL,
    -1,
    module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit_pairs(void) {
    import_array();
    return PyModule_Create(&moduledef);
}
