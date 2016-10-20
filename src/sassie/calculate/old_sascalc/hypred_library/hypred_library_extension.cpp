#include <math.h>
#include <Python.h>
#include <stdio.h>
#include <iostream>
#include "numpy/arrayobject.h"

#include "Atom.h"
#include "coor.h"
#include "box.h"
#include "pRDF.h"
#include "util.h"


// construct the prdf profile
static PyObject *func_construct_prdf(PyObject *self, PyObject *args)
{
	// local variables
	double cube_length;
	int nx, ny, nz;
	double x0, y0, z0;
	double prdf_spacing;
	int prdf_n;
	double cutoff_distance;
	size_t i;
	PyObject * item;

	// parse the argument list
    PyObject * py_solute_atom_types, * py_solvent_atom_types;
	PyObject * py_solute_atom_surf_mask;
    PyObject * py_solute_atom_coors, * py_solvent_atom_coors;
    if (!PyArg_ParseTuple(args, "O!O!O!O!O!ddddiiidid", &PyList_Type, &py_solute_atom_types, &PyArray_Type, &py_solute_atom_coors, &PyList_Type, &py_solvent_atom_types, &PyArray_Type, &py_solvent_atom_coors, &PyArray_Type, &py_solute_atom_surf_mask, &x0, &y0, &z0, &cube_length, &nx, &ny, &nz, &prdf_spacing, &prdf_n, &cutoff_distance))
		return NULL;

	// define types
	typedef Atom::SoluteAtom::SoluteAtomType solute_atom_type_t;
	typedef Atom::SolventAtom::SolventAtomType solvent_atom_type_t;

	// get the number of frames, and solute/solvent atoms
    PyArrayObject * pyArray_solute_atom_coors = (PyArrayObject*)PyArray_ContiguousFromObject(py_solute_atom_coors,PyArray_DOUBLE, 2,2);
    PyArrayObject * pyArray_solvent_atom_coors = (PyArrayObject*)PyArray_ContiguousFromObject(py_solvent_atom_coors,PyArray_DOUBLE, 3,3);
	const size_t tot_frames=pyArray_solvent_atom_coors->dimensions[0];
	const size_t n_solute_atoms = PyList_Size(py_solute_atom_types);
	const size_t n_solvent_atoms = PyList_Size(py_solvent_atom_types);
std::cout<<__LINE__<<std::endl;

	// get the solute atom surface mask
    PyArrayObject * pyArray_solute_atom_surf_mask = (PyArrayObject*)PyArray_ContiguousFromObject(py_solute_atom_surf_mask,PyArray_INT, 1,1);
std::cout<<__LINE__<<std::endl;

	// convert python types to C types
	solute_atom_type_t * solute_atom_types = new solute_atom_type_t[n_solute_atoms];
	for (i=0; i<n_solute_atoms; i++)
	{
		item  = PyList_GetItem(py_solute_atom_types, i);
		if (!PyString_Check(item)) prdf::utils::Error("solute atom type should be a string!");
		solute_atom_types[i] = PyString_AsString(item);
	}
	solvent_atom_type_t * solvent_atom_types = new solvent_atom_type_t[n_solvent_atoms];
	for (i=0; i<n_solvent_atoms; i++)
	{
		item  = PyList_GetItem(py_solvent_atom_types, i);
		if (!PyString_Check(item)) prdf::utils::Error("solvent atom type should be a string!");
		solvent_atom_types[i] = PyString_AsString(item);
	}
	double * solute_atom_coors = (double*)(pyArray_solute_atom_coors->data);
	double * solvent_atom_coors = (double*)(pyArray_solvent_atom_coors->data);
	int * solute_atom_surf_mask = (int*)(pyArray_solute_atom_surf_mask->data);
std::cout<<__LINE__<<std::endl;
	/*
	std::cout<<"c0: "<<x0<<" "<<y0<<" "<<z0<<std::endl;
	std::cout<<"nx/y/z: "<<nx<<" "<<ny<<" "<<nz<<std::endl;
	std::cout<<"cube_l: "<<cube_length<<std::endl;
	std::cout<<"prfd_n/spacing: "<<prdf_n<<" "<<prdf_spacing<<std::endl;
	std::cout<<"n_slute: "<<n_solute_atoms<<" n_solv: "<<n_solvent_atoms<<std::endl;
	std::cout<<"c0 solute: "<<solute_atom_coors[0]<<" "<<solute_atom_coors[1]<<" "<<solute_atom_coors[2]<<std::endl;
	std::cout<<"c1 solute: "<<solute_atom_coors[3]<<" "<<solute_atom_coors[4]<<" "<<solute_atom_coors[5]<<std::endl;
	std::cout<<"c0 solvent: "<<solvent_atom_coors[0]<<" "<<solvent_atom_coors[1]<<" "<<solvent_atom_coors[2]<<std::endl;
	std::cout<<"c1 solvent: "<<solvent_atom_coors[3]<<" "<<solvent_atom_coors[4]<<" "<<solvent_atom_coors[5]<<std::endl;
	for (i=0; i<n_solute_atoms; i++) std::cout<<solute_atom_surf_mask[i]<<std::endl;
	*/

	// box parameters
	Coor coor_box_origin(x0, y0, z0);
	// set up the box
	Box box(coor_box_origin, nx, ny, nz, cube_length);
std::cout<<__LINE__<<std::endl;

	// prdf parameters
	// set up prdf
	pRDF rdfs(prdf_n, prdf_spacing);
std::cout<<__LINE__<<std::endl;

	// set up the cubes in box
	box.setupCube(n_solute_atoms, solute_atom_types, solute_atom_coors, solute_atom_surf_mask, cutoff_distance); // set up the _p_cube_type_distance from the passed solute atom set
	// print out the types and distances values of all the not-null cubes for debugging
	box.fprintfTypeDistancePDB("box_type_distance.pdb");
std::cout<<__LINE__<<std::endl;

	// accumulate the solvent density into cubes
	for (size_t frame = 0; frame<tot_frames; ++frame)
	{
		// accumulate the solvent density into cubes
		box.accumulateCubeDensity(n_solvent_atoms, solvent_atom_types, solvent_atom_coors+frame*n_solvent_atoms*3);
	}
std::cout<<__LINE__<<std::endl;

	// finally update the pRDF based on the box
	rdfs.update(box);
std::cout<<__LINE__<<std::endl;

	// output pRDF
	std::ofstream fout("pRDF.txt");
	rdfs.print(fout,box,tot_frames);
	fout.close();
std::cout<<__LINE__<<std::endl;

	// clean up
	delete [] solute_atom_types;
	delete [] solvent_atom_types;
std::cout<<__LINE__<<std::endl;

	// return
	Py_INCREF(Py_None);
std::cout<<__LINE__<<std::endl;
	return Py_None;
}

/// predict the prdf profile
static PyObject *func_predict_prdf(PyObject *self, PyObject *args)
{
	// local variables
	double cube_length;
	int nx, ny, nz;
	double x0, y0, z0;
	double cutoff_distance;
	size_t i;
	PyObject * item;

	// parse the argument list
    PyObject * py_solute_atom_types;
	PyObject * py_solute_atom_surf_mask;
    PyObject * py_solute_atom_coors;
    if (!PyArg_ParseTuple(args, "O!O!O!ddddiiid", &PyList_Type, &py_solute_atom_types, &PyArray_Type, &py_solute_atom_coors, &PyArray_Type, &py_solute_atom_surf_mask, &x0, &y0, &z0, &cube_length, &nx, &ny, &nz, &cutoff_distance))
		return NULL;

	// define types
	typedef Atom::SoluteAtom::SoluteAtomType solute_atom_type_t;
	typedef Atom::SolventAtom::SolventAtomType solvent_atom_type_t;

	// get the number of frames, and solute/solvent atoms
    PyArrayObject * pyArray_solute_atom_coors = (PyArrayObject*)PyArray_ContiguousFromObject(py_solute_atom_coors,PyArray_DOUBLE, 2,2);
	const size_t n_solute_atoms = PyList_Size(py_solute_atom_types);

	// get the solute atom surface mask
    PyArrayObject * pyArray_solute_atom_surf_mask = (PyArrayObject*)PyArray_ContiguousFromObject(py_solute_atom_surf_mask,PyArray_INT, 1,1);

	// convert python types to C types
	solute_atom_type_t * solute_atom_types = new solute_atom_type_t[n_solute_atoms];
	for (i=0; i<n_solute_atoms; i++)
	{
		item  = PyList_GetItem(py_solute_atom_types, i);
		if (!PyString_Check(item)) prdf::utils::Error("solute atom type should be a string!");
		solute_atom_types[i] = PyString_AsString(item);
	}
	double * solute_atom_coors = (double*)(pyArray_solute_atom_coors->data);
	int * solute_atom_surf_mask = (int*)(pyArray_solute_atom_surf_mask->data);

	// box parameters
	Coor coor_box_origin(x0, y0, z0);
	// set up the box
	Box box(coor_box_origin, nx, ny, nz, cube_length);

	// instantiate pRDF from the external file
	pRDF rdfs("pRDF.txt");

	// set up the cubes in box
	box.setupCube(n_solute_atoms, solute_atom_types, solute_atom_coors, solute_atom_surf_mask, cutoff_distance); // set up the _p_cube_type_distance from the passed solute atom set
	// print out the types and distances values of all the not-null cubes for debugging
	box.fprintfTypeDistancePDB("box_type_distance.pdb");

	// predict the solvent density in cubes
	box.predictCubeDensity(rdfs, cutoff_distance);

	// output box
	box.fprintfTypeDensityPDB("box_type_density.pdb");

	// clean up
	delete [] solute_atom_types;

	// return
	Py_INCREF(Py_None);
	return Py_None;
}


static PyMethodDef Methods[] = 
{
	{ "construct_prdf", func_construct_prdf, METH_VARARGS },
	{ "predict_prdf", func_predict_prdf, METH_VARARGS },
	{ NULL, NULL }
} ;

PyMODINIT_FUNC initprdf()
{
    import_array();
	Py_InitModule("prdf", Methods);
}
