#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>

static PyObject* overlap(PyObject* self, PyObject* args) {
    PyArrayObject *coor1_array, *cutoffarray_array;
    double fact, x1, y1, z1, x2, y2, z2, diff2, dist;
    int check = 0, count = 1;
    npy_intp natoms1, npairs;

    if (!PyArg_ParseTuple(args, "O!O!d", &PyArray_Type, &coor1_array, &PyArray_Type, &cutoffarray_array, &fact)) {
        return NULL;
    }

    natoms1 = PyArray_DIM(coor1_array, 0);
    npairs = PyArray_DIM(cutoffarray_array, 0);

    for (npy_intp i = 0; i < natoms1 - 1; i++) {
        x1 = *(double *)PyArray_GETPTR2(coor1_array, i, 0);
        y1 = *(double *)PyArray_GETPTR2(coor1_array, i, 1);
        z1 = *(double *)PyArray_GETPTR2(coor1_array, i, 2);
        for (npy_intp j = i + 1; j < natoms1; j++) {
            x2 = *(double *)PyArray_GETPTR2(coor1_array, j, 0);
            y2 = *(double *)PyArray_GETPTR2(coor1_array, j, 1);
            z2 = *(double *)PyArray_GETPTR2(coor1_array, j, 2);
            diff2 = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2);
            dist = sqrt(diff2);
            if (dist < fact * *(double *)PyArray_GETPTR1(cutoffarray_array, count)) {
                check = 1;
                break;
            }
            count++;
        }
        if (check == 1) {
            break;
        }
    }

    return Py_BuildValue("i", check);
}

static PyMethodDef module_methods[] = {
    {"overlap", overlap, METH_VARARGS, "Calculate overlap"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "vdw_overlap",
    NULL,
    -1,
    module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit_vdw_overlap(void) {
    import_array();
    return PyModule_Create(&moduledef);
}
