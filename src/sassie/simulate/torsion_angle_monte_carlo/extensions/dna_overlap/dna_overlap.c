#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>

static PyObject* wca_d(PyObject* self, PyObject* args) {
    PyArrayObject *coor, *wca;
    int tbead, n;
    double w;
    if (!PyArg_ParseTuple(args, "O!iidO!", &PyArray_Type, &coor, &tbead, &w, &n, &PyArray_Type, &wca)) {
        return NULL;
    }

    double cutoff = pow(2.0, 1.0/6.0) * w;
    double *coor_data = (double*) PyArray_DATA(coor);
    double *wca_data = (double*) PyArray_DATA(wca);

    for (int i = tbead; i < n; i++) {
        double x1 = coor_data[i*3];
        double y1 = coor_data[i*3 + 1];
        double z1 = coor_data[i*3 + 2];
        for (int j = 0; j < i; j++) {
            double x2 = coor_data[j*3];
            double y2 = coor_data[j*3 + 1];
            double z2 = coor_data[j*3 + 2];

            double dx = x2 - x1;
            double dy = y2 - y1;
            double dz = z2 - z1;
            double rij = sqrt(dx*dx + dy*dy + dz*dz);

            if (rij < cutoff) {
                wca_data[i*n + j] = pow(w/rij, 12) - pow(w/rij, 6) + 0.25;
                wca_data[j*n + i] = wca_data[i*n + j]; // Symmetric matrix
            } else {
                wca_data[i*n + j] = 0.0;
                wca_data[j*n + i] = 0.0;
            }
        }
    }

    Py_RETURN_NONE;
}

static PyObject* overlap_skip_n(PyObject* self, PyObject* args) {
    PyArrayObject *coor;
    int natoms, nskip;
    double cutoff;
    if (!PyArg_ParseTuple(args, "O!iid", &PyArray_Type, &coor, &natoms, &nskip, &cutoff)) {
        return NULL;
    }

    double *coor_data = (double*) PyArray_DATA(coor);
    int check = 0;
    int i_index = 0, j_index = 0;

    for (int i = 0; i < natoms - nskip - 1; i++) {
        double x1 = coor_data[i*3];
        double y1 = coor_data[i*3 + 1];
        double z1 = coor_data[i*3 + 2];
        for (int j = i + nskip + 1; j < natoms; j++) {
            double x2 = coor_data[j*3];
            double y2 = coor_data[j*3 + 1];
            double z2 = coor_data[j*3 + 2];

            double diff2 = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2);
            double dist = sqrt(diff2);

            if (dist < cutoff) {
                check = 1;
                i_index = i + 1; // Fortran is 1-indexed
                j_index = j + 1; // Adjust for Python's 0-indexing
                goto endloop;
            }
        }
    }
endloop:

    return Py_BuildValue("iii", check, i_index, j_index);
}

static PyObject* overlap2(PyObject* self, PyObject* args) {
    PyArrayObject *coor1a, *coor1b;
    double cutoff;
    if (!PyArg_ParseTuple(args, "O!O!d", &PyArray_Type, &coor1a, &PyArray_Type, &coor1b, &cutoff)) {
        return NULL;
    }

    int n1a = PyArray_DIM(coor1a, 0);
    int n1b = PyArray_DIM(coor1b, 0);
    double *coor1a_data = (double*) PyArray_DATA(coor1a);
    double *coor1b_data = (double*) PyArray_DATA(coor1b);
    int check = 0;
    int i_index = 0, j_index = 0;

    for (int i = 0; i < n1a; i++) {
        double x1 = coor1a_data[i*3];
        double y1 = coor1a_data[i*3 + 1];
        double z1 = coor1a_data[i*3 + 2];
        for (int j = 0; j < n1b; j++) {
            double x2 = coor1b_data[j*3];
            double y2 = coor1b_data[j*3 + 1];
            double z2 = coor1b_data[j*3 + 2];

            double diff2 = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2);
            double dist = sqrt(diff2);

            if (dist < cutoff) {
                check = 1;
                i_index = i + 1; // Adjust for Python's 0-indexing to match Fortran's 1-indexing
                j_index = j + 1;
                goto endloop;
            }
        }
    }
endloop:

    return Py_BuildValue("iii", check, i_index, j_index);
}

static PyMethodDef module_methods[] = {
    {"overlap2", (PyCFunction)overlap2, METH_VARARGS, "Check for overlap between two sets of coordinates within a cutoff distance"},
    {"wca_d", (PyCFunction)wca_d, METH_VARARGS, "Calculate interaction energy for a given set of coordinates"},
    {"overlap_skip_n", (PyCFunction)overlap_skip_n, METH_VARARGS, "Check for overlap with skip within a cutoff distance"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "dna_overlap", /* name of module */
    NULL,              /* module documentation, may be NULL */
    -1,                /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    module_methods
};

PyMODINIT_FUNC PyInit_dna_overlap(void) {
    import_array(); // Necessary for NumPy initialization
    return PyModule_Create(&moduledef);
}