import numpy
# import sassie.simulate.monte_carlo.monte_carlo_utilities.double_stranded_nucleic.create_sasmol as create_sasmol
import create_sasmol

def create_sasmol_as_cpp(mol):

    index = numpy.array(mol.index())
    name = numpy.array(mol.name(),dtype = 'string').tolist()
    loc = numpy.array(mol.loc(),dtype = 'string').tolist()
    resname = numpy.array(mol.resname(),dtype = 'string').tolist()
    chain = numpy.array(mol.chain(),dtype = 'string').tolist()
    resid = numpy.array(mol.resid())
    rescode = numpy.array(mol.rescode(),dtype = 'string').tolist()
    x = numpy.array(mol.coor()[:,:,0]).flatten().copy()
    y = numpy.array(mol.coor()[:,:,1]).flatten().copy()
    z = numpy.array(mol.coor()[:,:,2]).flatten().copy()
    occupancy = numpy.array(mol.occupancy(),dtype = 'string').tolist()
    beta = numpy.array(mol.beta(),dtype = 'string').tolist()
    segname = numpy.array(mol.segname(),dtype = 'string').tolist()
    element = numpy.array(mol.element(),dtype = 'string').tolist()
    charge = numpy.array(mol.charge(),dtype = 'string').tolist()

    cpp_object = create_sasmol.create_and_return_cpp_class_instance(mol.natoms(), mol.number_of_frames(), index, name, loc, resname, chain, resid, rescode, x, y, z, occupancy, beta, segname, element, charge)
    return cpp_object
