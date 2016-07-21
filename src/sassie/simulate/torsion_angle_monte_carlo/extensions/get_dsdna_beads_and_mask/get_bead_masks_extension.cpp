#include "Python.h"
#include "numpy/arrayobject.h"

#include "sasmol.h"
#include "sassubset.h"

#include <iostream>
#include <sstream>


///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
PyObject *get_bead_masks_in_cpp(PyObject *self, PyObject *args){

	PyObject *python_sasmol= NULL;
	PyObject *python_bp_in_each_bead= NULL;
    int number_beads;
    int chain1_res_a;
    int chain2_res_b;
    int number_base_pairs;
    PyObject * python_nucleic_segname1;
    PyObject * python_nucleic_segname2;

    int chain1_res_b;
    int chain2_res_a;
	
	if (!PyArg_ParseTuple(args, "OOiiiiOO", &python_sasmol, &python_bp_in_each_bead, &number_beads, &chain1_res_a, &chain2_res_b, &number_base_pairs, &python_nucleic_segname1, &python_nucleic_segname2))
		return NULL;

    sasmol::SasMol * p = static_cast<sasmol::SasMol*>(PyCObject_AsVoidPtr(python_sasmol));
    int natoms = p->_natoms();

    PyArrayObject * pyArray_bp_in_each_bead = (PyArrayObject*)PyArray_ContiguousFromObject(python_bp_in_each_bead,PyArray_INT64, 1,1);
    long long *bp_in_each_bead = (long long*)(pyArray_bp_in_each_bead->data);

    std::string nucleic_segname1 = PyString_AsString(python_nucleic_segname1);
    std::string nucleic_segname2 = PyString_AsString(python_nucleic_segname2);

    std::string bead_filter;
    std::stringstream ss_chain1_bead_filter, ss_chain2_bead_filter;

    int dims[2]; dims[0]=number_beads; dims[1]=natoms;
    PyArrayObject * bead_masks = (PyArrayObject*)PyArray_FromDims(2, dims, PyArray_INT32);
    char * p_data;

    for (int j=0; j<number_beads; ++j)
    {
        if (bp_in_each_bead[j] > 1)
        {
            // create the basis filters fore ach bead
            if (j+1 == number_beads)
            {
                // accomodate for a non-divisible number of residues
                bp_in_each_bead[j] = number_base_pairs - bp_in_each_bead[j] * j;
            }
                
            // Get the atoms from DNA strand 1
            chain1_res_b = chain1_res_a + bp_in_each_bead[j];
            ss_chain1_bead_filter.str("");
            ss_chain1_bead_filter<<"((resid "<<chain1_res_a<<" "<<chain1_res_b-1<<") and (tag "<<nucleic_segname1<<")) or ";
            chain2_res_a = chain2_res_b - bp_in_each_bead[j];
            ss_chain2_bead_filter.str("");
            ss_chain2_bead_filter<<"((resid "<<chain2_res_a<<" " + chain2_res_b<<") and (tag "<< nucleic_segname2<<"))";

            // setup for next iteration
            chain1_res_a = chain1_res_b;
            chain2_res_b = chain2_res_a;
        }
        else
        {
            ss_chain1_bead_filter.str("");
            ss_chain1_bead_filter<<"((resid "<<chain1_res_a<<") and (tag "<<nucleic_segname1<<")) or ";
            ss_chain2_bead_filter.str("");
            ss_chain2_bead_filter<<"((resid "<<chain2_res_b<<") and (tag "<<nucleic_segname2<<"))";

            // setup for next iteration
            chain1_res_a += 1;
            chain2_res_b -= 1;
        }
            
        bead_filter = ss_chain1_bead_filter.str()+ss_chain2_bead_filter.str();
        
        // create the bead masks to select the atoms from the aa mol
        //std::cout<<"bead_filter: "<<bead_filter<<std::endl;
        boost::dynamic_bitset<> mask = p->get_subset_mask(bead_filter);
        //std::cout<<mask<<std::endl;
        //handle_error(other_self, error)
        
        // store the mask for the reverse coarse-graining
        for (int k=0; k<natoms; ++k)
        {
            p_data = bead_masks->data + j*bead_masks->strides[0] + k*bead_masks->strides[1];
            *p_data = int(mask[k]);
        }
    }

    return PyArray_Return(bead_masks);
}

///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
static PyMethodDef exampleMethods[] = {
	{ "get_bead_masks_in_cpp", get_bead_masks_in_cpp, METH_VARARGS },
	{ NULL, NULL }
} ;


PyMODINIT_FUNC
initget_bead_masks(){
	PyObject *m, *m1;
    m = Py_InitModule("get_bead_masks", exampleMethods);
	m1 = Py_InitModule("get_bead_masks_in_cpp", exampleMethods);
	import_array();
}

