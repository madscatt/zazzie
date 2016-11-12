import sassie.calculate.sascalc.sascalc_library.sascalc_api as sascalc_api
import box_converge_lib
import numpy, glob, os

class sascalc_results:
    pass

class SasCalc:
    def __init__(self):
        self.current_frame = 0

    def initialize(self,sascalc_inputs):
        self.sascalc_inputs = sascalc_inputs
        sascalc_inputs.coor = numpy.transpose(sascalc_inputs.coor,axes=(0,2,1))
        self.o = sascalc_api.initialize(sascalc_inputs)
        return self

    def batch_load(self):
        inputs = self.sascalc_inputs
        frames_per_batch = inputs.frames_per_batch
        number_of_frames = inputs.number_of_frames
        if self.current_frame == 0:
            offset = 0
        else:
            offset = frames_per_batch
        extend = min(frames_per_batch, number_of_frames-self.current_frame)
        dummy = sascalc_api.batch_load(self.o,offset,extend)
        self.current_frame += extend 

    def calculate(self,frame):
        inputs = self.sascalc_inputs
        results = sascalc_results()
        results.Iq_neutron_array = numpy.zeros((3*3, inputs.B_neutron_array.shape[1], inputs.Q.size)) # vacuum, solvent, complete,  vacuum_real, solvent_real, complete_real, vacuum_imag, solvent_imag, complete_imag
        results.Iq_xray_array = numpy.zeros((3*3, inputs.B_xray_array.shape[1], inputs.Q.size)) # vacuum, solvent, complete,  vacuum_real, solvent_real, complete_real, vacuum_imag, solvent_imag, complete_imag
        results.Ngv_converged_neutron_array = numpy.zeros((3, inputs.B_neutron_array.shape[1]),dtype=numpy.int32)
        results.Ngv_converged_xray_array = numpy.zeros((3, inputs.B_xray_array.shape[1]),dtype=numpy.int32)
        results.Pr= numpy.zeros(inputs.R.size)
        dummy = sascalc_api.calculate(self.o,frame,inputs,results)
        return results

    def box_converge(self):
        inputs = self.sascalc_inputs
        output_folder = inputs.output_folder
        box_converge_lib.converge_real_space(inputs.mol, output_folder)
        if inputs.xon in ['neutron','neutron_and_xray']:
            output_dir = glob.glob(os.path.join(output_folder,'neutron_D2Op_[0-9]*'))[0]
        else:
            output_dir = os.path.join(output_folder,'xray')
        box_converge_lib.converge_sas_space(inputs.runname, output_dir)
    
    def clean(self):
        sascalc_api.clean(self.o)
