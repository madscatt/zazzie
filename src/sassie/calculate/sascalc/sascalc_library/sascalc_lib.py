import sassie.calculate.sascalc.sascalc_library.sascalc_api as sascalc_api
import numpy, glob, os

class ScResults:
    
    def __init__(self, results, runname):
        '''
        initiliazation
        '''
        self.results = results
        self.runname = runname

    def save(self, name):
        '''
        save results
        '''
        sascalc_api.save(self.results, name, self.runname)

class SasCalc:

    def __init__(self, mol, mvars, scvars):
        '''
        initialization
        '''
        if mvars.xon in ['neutron','neutron_and_xray']:
            B_neutron_array = scvars.B_neutron_array[:,0,:]
        else:
            B_neutron_array = numpy.array([])
        if mvars.xon in ['xray','neutron_and_xray']:
            B_xray_array = scvars.B_xray_array[:,0,:,:]
        else:
            B_xray_array = numpy.array([])

        self.runname = scvars.runname
        self.o = sascalc_api.initialize(mol.natoms(), mvars, B_neutron_array, B_xray_array)

    def calculate(self, coordinates, Nframes):
        '''
        calculator
        '''
        results = sascalc_api.calculate(self.o, coordinates, Nframes)
        return ScResults(results, self.runname)

    def epilogue(self):
        sascalc_api.clean(self.o)

