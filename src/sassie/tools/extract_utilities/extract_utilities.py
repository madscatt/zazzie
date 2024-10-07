'''
    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

#
#	EXTRACT UTILITIES
#
#	08/27/2012	--	initial coding			                        :	jc
#	09/06/2012	--	fixed for large DCD files	                    : 	hz
#	04/19/2015	--	added SAS options                               :   jc
#	10/09/2016	--	added sascalc and existing folder timestamps    :   jc
#
#	 1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#								       *      **

import os
import sys
import locale
import numpy
import time
import glob
import sassie.interface.input_filter as input_filter
import sasmol.system as system
import sassie.util.module_utilities as module_utilities
import sassie.util.sasconfig as sasconfig
import sassie.util.folder_management as folder_management
#import folder_management as folder_management

import datetime

'''
	EXTRACT_UTILITIES is a module that allows one to extract coordinates and/or
    SAS profiles from PDB/DCD files and/or a directory containing SAS profiles.
	The multi-frame files trajectory is saved into new a PDB/DCD file.  SAS profiles
    are saved to a new directory.   Options are provided to extract
	single structures and/or SAS, structures and/or SAS over a range, by reading a list of frame 
	numbers from a simple text file, via a 'weights' file, or by a user-supplied
    sampling frequency.


###
###	NOTE: should be straightforward to add basis filtering logic to pull out
###	      subsets.  Although this will take some work and may be better to
###           include with "builder" and/or "pdbrx" type projects
###
###	NOTE2: You could re-write this so that the single frame and range options
###	       work just like the text_file and weight_file versions (i.e. using masks)
###
###
'''

if sasconfig.__level__ == "DEBUG":
    DEBUG = True

app = 'extract_utilities'

class module_variables():

    def __init__(self, parent=None):
        self.app = app

class extract_utilities_input_variables():

    def __init__(self, parent=None):
        pass


class extract_utilities():

    def __init__(self, parent=None):
        pass

    def main(self, input_variables, txtOutput):

        self.module_variables = module_variables()

        self.extract_utilities_input_variables = extract_utilities_input_variables()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.run_utils.general_setup(self)

        self.initialization()

        self.extract_utilities()

        self.epilogue()

        return

#    pgui performs this function
#    def print_failure(message, txtOutput):
#
#        txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
#        txtOutput.put(">>>> RUN FAILURE <<<<\n")
#        txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
#        txtOutput.put(message)
#
#        return


    def unpack_variables(self,variables):
        '''
        method to extract variables into system wise class instance
        '''

        log = self.log
        mvars = self.module_variables 
        log.debug('in unpack_variables')
        
        mvars.run_name = variables['run_name'][0]
        mvars.path = variables['path'][0]
        mvars.pdb_filename = variables['pdb_filename'][0]
        mvars.trajectory_filename = variables['trajectory_filename'][0]
        mvars.option = variables['option'][0]
        mvars.local_value = variables['local_value'][0]
        mvars.output_filename = variables['output_filename'][0]
        mvars.extract_trajectory = variables['extract_trajectory'][0]
        mvars.extract_sas = variables['extract_sas'][0]
        mvars.sas_type = variables['sas_type'][0]
        mvars.sas_paths = variables['sas_paths'][0]

        log.debug(vars(mvars))

        return 


    def copy_spec_files(self,sas_files,suffix):

        log = self.log
        mvars = self.module_variables
        avars = self.extract_utilities_input_variables

        string_fill = 5
        number = 1
        number_of_saved_files = 0

        for name in sas_files:
            if str(number) in avars.mask:
                number_of_saved_files += 1
                numst = str(number_of_saved_files).zfill(string_fill)
                runst = mvars.run_name + '_' + numst + suffix[1:]
                log.debug('runst: '+ runst)
                this_file = os.path.join(avars.sas_output_path, runst)
                cpst = 'cp  ' + name + ' ' + this_file
                log.debug('cpst: '+ cpst)
                os.system(cpst)

            number += 1

        return number_of_saved_files


    def get_range_sas(self):

        mvars = self.module_variables

        local_values = mvars.local_value.split('-')
        first = locale.atoi(local_values[0])
        last = locale.atoi(local_values[1])

        return list(range(first, last + 1))


    def get_text_file_sas(self):

        mvars = self.module_variables

        infile = open(mvars.local_value, 'r').readlines()
        mask = []
        for i in range(len(infile)):
            lin = infile[i].split()
            if(len(lin) > 0):
                mask.append(lin[0])

        return mask


    def get_weight_file_sas(self):

        log = self.log
        mvars = self.module_variables

        mask = []
        weights = numpy.loadtxt(mvars.local_value)
        fweights = weights[:, 2]
        file_numbers = weights[:, 0]
        x2 = weights[:, 1]
        for i in range(len(x2)):
            if(fweights[i] == 1):
                log.debug('st = '+str(weights[i][0])+' : x2 = '+str(x2[i]))
                pass
        frame_list = ((numpy.nonzero(fweights))[0]).tolist()
        mask = [str(int(file_numbers[i])) for i in frame_list]

        return mask


    def get_frequency(self, number_of_spec_files, **kwargs):

        log = self.log
        mvars = self.module_variables

        mask = []
        log.debug('nspecfil: '+ str(number_of_spec_files))

        try:
            coordinate_flag = kwargs['coordinate_flag']
            print('coordinate flag: ',coordinate_flag)
            print('extracting coordinates: local_value = ', mvars.local_value)
            coordinate_flag = True

        except:
            coordinate_flag = False

        print('coordinate flag after check: ', coordinate_flag)

        if not coordinate_flag:
            for number in range(1, number_of_spec_files + 1, int(mvars.local_value)):
                mask.append(str(number))
        elif coordinate_flag:
            for number in range(0, number_of_spec_files, int(mvars.local_value)):
                mask.append(number)

        return mask


    def get_sas_mask(self,**kwargs):

        log = self.log
        mvars = self.module_variables
        log.debug('getting SAS mask')

        mask = []

        if mvars.option == 'single_frame':
            mask.append(mvars.local_value)
        elif mvars.option == 'range':
            rangelist = self.get_range_sas()
            mask = [str(i) for i in rangelist]
        elif mvars.option == 'text_file':
            mask = self.get_text_file_sas()
        elif mvars.option == 'weight_file':
            mask = self.get_weight_file_sas()
        elif mvars.option == 'sampling_frequency':
            number_of_spec_files = kwargs['number_of_spec_files']
            mask = self.get_frequency(number_of_spec_files)
        elif mvars.option == 'all':
            number_of_spec_files = kwargs['number_of_spec_files']
            mvars.local_value = '1'
            mask = self.get_frequency(number_of_spec_files)

        log.debug('mask: '+str(mask))

        return mask

    def get_single_frame(self):

        mvars = self.module_variables
    
        frame = locale.atoi(mvars.local_value) - 1
        return [frame]


    def get_range(self):

        mvars = self.module_variables

        local_values = mvars.local_value.split('-')
        first = locale.atoi(local_values[0]) - 1
        last = locale.atoi(local_values[1]) - 1
        return list(range(first, last + 1))


    def get_text_file(self):

        log = self.log
        mvars = self.module_variables
        
        infile = open(mvars.local_value, 'r').readlines()

        frame_list = []
        for i in range(len(infile)):
            lin = infile[i].split()
            if(len(lin) > 0):
                this_value = locale.atoi(lin[0]) - 1
                frame_list.append(this_value)
        log.debug('frame_list: '+str(frame_list))

        return frame_list


    def get_weight_file(self):

        log = self.log
        mvars = self.module_variables
    
        weights = numpy.loadtxt(mvars.local_value)
        fweights = weights[:, 2]
        x2 = weights[:, 1]
        for i in range(len(x2)):
            if(fweights[i] == 1):
                log.debug('st = '+str(weights[i][0])+' : x2 = '+str(x2[i]))
                pass
        frame_list = ((numpy.nonzero(fweights))[0]).tolist()
        log.debug('frame_list[0]: '+str(frame_list[0]))
        return frame_list


    def get_number_of_frames(self):

        mvars = self.module_variables
        avars = self.extract_utilities_input_variables

        m = system.Molecule(0)

        if avars.infile_type == 'pdb':
            m.read_pdb(avars.trajectory_filename)
            number_of_frames = m.number_of_frames()
        elif avars.infile_type == 'dcd':
            dcdfile = m.open_dcd_read(avars.trajectory_filename)
            number_of_frames = dcdfile[2]

        return number_of_frames



    def initialization(self):
        '''
        method to prepare for extraction of trajectories and/or SAS files
        '''

        log = self.log
        log.debug('in initialization')
        pgui = self.run_utils.print_gui

        mvars = self.module_variables
        avars = extract_utilities_input_variables

        log.debug(vars(mvars))
        log.debug(vars(avars))

        ''' directory and file preparation '''

        avars.pdb_filename = mvars.path + mvars.pdb_filename
        log.debug('pdb filename: '+avars.pdb_filename)
        avars.trajectory_filename = mvars.path + mvars.trajectory_filename
        log.debug('trajectory filename: '+avars.trajectory_filename)

        avars.output_path = mvars.run_name + '/extract_utilities/'
        log.debug('output_path: '+avars.output_path)

#        vers = 'version 0.2 : 04/23/15 : jc'
        direxist = os.path.exists(avars.output_path)
        if(direxist == 0):
            try:
                result = os.system('mkdir -p ' + avars.output_path)
            except:
                message = 'can not create project directory: ' + avars.output_path
                message += '\nstopping here\n'
                pgui(message)
                sys.exit(1)
            if(result != 0):
                message = 'can not create project directory: ' + avars.output_path
                message += '\nstopping here\n'
                pgui(message)
                sys.exit(1)


        if mvars.extract_trajectory:

            avars.infile_type = avars.trajectory_filename[-3:]

            if(mvars.option == 'single_frame'):
                avars.rangelist = self.get_single_frame()

            elif(mvars.option == 'range'):
               avars.rangelist = self.get_range()

            elif(mvars.option == 'text_file'):
                avars.rangelist = self.get_text_file()

            elif(mvars.option == 'weight_file'):
                avars.rangelist = self.get_weight_file()

            elif(mvars.option == 'sampling_frequency'):
                number_of_frames = self.get_number_of_frames()
                avars.rangelist = self.get_frequency(number_of_frames, coordinate_flag=True)

            elif(mvars.option == 'all'):
                mvars.local_value = 1
                number_of_frames = self.get_number_of_frames()
                avars.rangelist = self.get_frequency(number_of_frames, coordinate_flag=True)

        if mvars.extract_sas:

            avars.sas_paths = [x.strip() for x in mvars.sas_paths.split(',')]
            base_paths = []
            for sas_path in avars.sas_paths:
                pgui("sas_path = %s\n" % (sas_path))
                base_paths.append(os.path.basename(os.path.normpath(sas_path)))

            avars.sas_output_paths = []

            if(mvars.sas_type == 0):
                for base_path in base_paths:
                    avars.sas_output_paths.append(os.path.join(
                        avars.output_path, 'sascalc', base_path))
                avars.suffix = ['*.iq', '*.log']
                avars.extra = []
                avars.num_ex_files = [0, 0, 0, 0]

            elif(mvars.sas_type == 1):
                avars.sas_output_paths.append(os.path.join(avars.output_path, 'xtal2sas'))
                avars.suffix = ['*.iq', '*.log']
                avars.extra = ['*.inf', '*.crd', '*.ans', '*.pr']
                avars.num_ex_files = [1, 1, 1, 1]
            elif(mvars.sas_type == 2):
                avars.sas_output_paths.append(os.path.join(avars.output_path, 'cryson'))
                avars.suffix = ['*.int', '*.log']
                avars.extra = ['*.sav', '*.flm', '*.alm', '*.ans']
                avars.num_ex_files = [1, 1, 1, 1]
            elif(mvars.sas_type == 3):
                avars.sas_output_paths.append(os.path.join(avars.output_path, 'crysol'))
                avars.suffix = ['*.int', '*.log']
                avars.extra = ['*.sav', '*.flm', '*.alm', '*.ans']
                avars.num_ex_files = [1, 1, 1, 1]
   
               
        return

    def extract_utilities(self):

        '''
        EXTRACT_UTILITIES is a module that allows one to extract coordinates and/or
        SAS profiles from PDB/DCD files and/or a directory containing SAS profiles.

	    INPUT:  variable descriptions
	
        run_name:                string      project name
        path:                   string      input/output filepath      
        pdb_filename            string      input pdb file
        trajectory_filename     string      input pdb or dcd file                          
        option:                 string      extract option (single_frame, range, all, weight_file, text_file)                 
        local_value:            string      value of option (frame value, range, all, weight_file name, text_file name)                
        output_filename:        string      output pdb or dcd file   
        extract_trajectory:     boolean     extract frames from trajectory (True or False)                
        extract_sas :           boolean     extract corresponding SAS files (True or False)                 
        sas_type:               integer     integer depends on SAS file type (SasCalc, Xtal2sas, Cryson, Crysol)    
        sas_path:               string      path to SAS files 


	    OUTPUT: files stored in "run_name"/extract_utilities directory:
        	
        PDB or DCD file with requested coordinates
            and/or
        Folder containing SAS profiles

        '''

        log = self.log
        pgui = self.run_utils.print_gui
        log.debug('in extract_utilities')

        mvars = self.module_variables
        avars = self.extract_utilities_input_variables

        log.debug(vars(mvars))
        log.debug(vars(avars))

        # ttxt=time.ctime()
        ttxt = time.asctime(time.gmtime(time.time()))
        st = ''.join(['=' for x in range(60)])

        pgui("\n%s \n" % (st))
        pgui("DATA FROM RUN: %s \n\n" % (ttxt))

        if mvars.extract_trajectory:

            pgui("reading frames from %s \n" % (avars.trajectory_filename))

            if mvars.extract_sas:
                fraction_done = 0.25
            else:
                fraction_done = 0.50

            report_string = 'STATUS\t' + str(fraction_done)
            pgui(report_string)

            self.extract_coords()

        if mvars.extract_sas:

            pgui("\nextracting SAS profiles\n")

            fraction_done = 0.5
            report_string = 'STATUS\t' + str(fraction_done)
            pgui(report_string)

            self.extract_sas_files()

        fraction_done = 1.0
        report_string = 'STATUS\t' + str(fraction_done)
        pgui(report_string)

        return


    def extract_coords(self):
        '''
        method to extract requested coordinates
        '''

        log = self.log
        pgui = self.run_utils.print_gui
        log.debug('extracting coords')

        mvars = self.module_variables
        avars = extract_utilities_input_variables

        log.debug(vars(mvars))
        log.debug(vars(avars))

        m1 = system.Molecule(0)
        m1.read_pdb(avars.pdb_filename)
        natoms = m1.natoms()

        output_filename = avars.output_path + mvars.output_filename

        coor = numpy.zeros((1, natoms, 3), numpy.float32)

        pgui("writing frames to %s \n" % (output_filename))
        if(avars.infile_type == 'dcd'):
            pgui('> input file is a dcd file')
            if(output_filename[-3:] == 'dcd'):
                pgui('> output file is a dcd file')
                j = 0
                m2 = system.Molecule(1)
                m2.read_pdb(avars.pdb_filename, fastread=True)
                coor = numpy.zeros((1, natoms, 3), numpy.float32)
                dcdfile = m1.open_dcd_read(avars.trajectory_filename)
                number_of_frames = dcdfile[2]
                pgui('number_of_frames: '+ str(number_of_frames))

                log.debug('rangelist: '+ str(avars.rangelist))
                log.debug('max(rangelist): '+ str(max(avars.rangelist)))
                print('rangelist: ', avars.rangelist)
                print('max(rangelist): ', max(avars.rangelist))
                dcdoutfile = m2.open_dcd_write(output_filename)
                for i in range(number_of_frames):
                    print('.', end=' ') 
                    sys.stdout.flush()                  
                    m1.read_dcd_step(dcdfile, i)
                    if i in avars.rangelist:
                        pgui('\nextracting coordinates from frame: '+str(i))
                        coor[0, :, :] = m1.coor()[0]
                        m2.setCoor(coor)
                        m2.write_dcd_step(dcdoutfile, 0, j + 1)
                        j += 1
                    if(i > max(avars.rangelist) + 1):
                        break

                m2.close_dcd_write(dcdoutfile)
                m1.close_dcd_read(dcdfile[0])

            elif(output_filename[-3:] == 'pdb'):
                pgui('> output file is a pdb file')
                m2 = system.Molecule(1)
                m2.read_pdb(avars.pdb_filename)  # ,fastread = True)
                j = 0
                dcdfile = m1.open_dcd_read(avars.trajectory_filename)
                number_of_frames = dcdfile[2]
                pgui('number_of_frames: '+ str(number_of_frames))
                for i in range(number_of_frames):
                    print('.', end=' ')
                    sys.stdout.flush()
                    m1.read_dcd_step(dcdfile, i)
                    coor[0, :, :] = m1.coor()[0]
                    m2.setCoor(coor)
                    if i in avars.rangelist:
                        pgui('\nextracting coordinates from frame: '+str(i))
                        if(j == 0):
                            m2.write_pdb(output_filename, 0, 'w', model=j + 1)
                        else:
                            m2.write_pdb(output_filename, 0, 'a', model=j + 1)
                        j += 1
                    if(i > max(avars.rangelist) + 1):
                        break

                with open(output_filename, "a") as myfile:
                    myfile.write("END\n")

        elif(avars.infile_type == 'pdb'):
            m1.read_pdb(avars.trajectory_filename)
            natoms = m1.natoms()
    #	    coor = numpy.zeros(((last-first+1),natoms,3),numpy.float32)
            coor = numpy.zeros((1, natoms, 3), numpy.float32)
            number_of_frames = m1.number_of_frames()

            pgui('> input file is a pdb file')
            if(output_filename[-3:] == 'dcd'):
                pgui('> output file is a dcd file')
                j = 0
                m2 = system.Molecule(1)
                m2.read_pdb(avars.trajectory_filename, fastread=True)
                log.debug('rangelist: '+ str(avars.rangelist))
                log.debug('max(rangelist): '+ str(max(avars.rangelist)))
                print('rangelist: ', avars.rangelist)
                print('max(rangelist): ', max(avars.rangelist))                
                dcdoutfile = m2.open_dcd_write(output_filename)
                for i in range(number_of_frames):
                    print('.', end=' ')
                    sys.stdout.flush()
                    if i in avars.rangelist:
                        pgui('\nextracting coordinates from frame: '+str(i))
                        coor[0, :, :] = m1.coor()[i]
                        m2.setCoor(coor)
                        m2.write_dcd_step(dcdoutfile, 0, j + 1)
                        j += 1
                    if(i > max(avars.rangelist) + 1):
                        break

            elif(output_filename[-3:] == 'pdb'):
                pgui('> output file is a pdb file')
                m2 = system.Molecule(1)
                m2.read_pdb(avars.pdb_filename, fastread=True)
                j = 0
                for i in range(number_of_frames):
                    print('.', end=' ')
                    sys.stdout.flush()
                    if i in avars.rangelist:
                        pgui('\nextracting coordinates from frame = '+str(i))

                        coor[0, :, :] = m1.coor()[i]
                        m2.setCoor(coor)
                        if(j == 0):
                            m2.write_pdb(output_filename, 0, 'w', model=j + 1)
                        else:
                            m2.write_pdb(output_filename, 0, 'a', model=j + 1)
                        j += 1
                    if(i > max(avars.rangelist) + 1):
                        break

                with open(output_filename, "a") as myfile:
                    myfile.write("END\n")

        pgui("wrote %i frames to %s \n" %
                    (len(avars.rangelist), output_filename))
        return


    def extract_sas_files(self):
        '''
        method to extract requested SAS files
        '''
        
        log = self.log
        pgui = self.run_utils.print_gui
        log.debug('extracting SAS')

        mvars = self.module_variables
        avars = extract_utilities_input_variables

        log.debug(vars(mvars))
        log.debug(vars(avars))

        for sas_output_path in avars.sas_output_paths:
            direxist = os.path.exists(sas_output_path)
            if(direxist):
                try:
                    modification_date = folder_management.modification_date(
                        sas_output_path)
                    new_directory_name = sas_output_path + '_' + modification_date + '_UTC'
                    result = os.system('mv ' + sas_output_path +
                                    ' ' + new_directory_name)
                    result = os.system('mkdir -p ' + sas_output_path)
                except:
                    message = 'can not create project directory: ' + sas_output_path
                    message = '\n or can not create backup project directory: ' + \
                        new_directory_name + '_UTC'
                    message += '\nstopping here\n'
                    pgui(message)
                    sys.exit()
            else:
                result = os.system('mkdir -p ' + sas_output_path)

        copy_extras = False

        num_iq_files = 1
        num_log_files = 1

        count = 0

        for sas_path in avars.sas_paths:
            avars.sas_output_path = avars.sas_output_paths[count]
            spec_files = []
            log_files = []
            extra_files = []

            for name in glob.glob(os.path.join(sas_path, avars.suffix[0])):
                spec_files.append(name)
            for name in glob.glob(os.path.join(sas_path, avars.suffix[1])):
                log_files.append(name)
            for ex in avars.extra:
                this_extra = []
                for name in glob.glob(os.path.join(sas_path, ex)):
                    this_extra.append(name)
                if(copy_extras == False and len(this_extra) > 0):
                    copy_extras = True
                    pgui('copying extra sas files')
                this_extra.sort()
                extra_files.append(this_extra)	
                

            spec_files.sort()
            log_files.sort()


            if(len(spec_files) == 0):
                log.debug('sas_path,suffix: ' + os.path.join(sas_path, avars.suffix[0]))
                sys.exit()

            if(mvars.option in ['sampling_frequency', 'all']):
                avars.mask = self.get_sas_mask(number_of_spec_files=len(spec_files))

            else:
                avars.mask = self.get_sas_mask()

            num_iq_files = self.copy_spec_files(spec_files,avars.suffix[0])
            num_log_files = self.copy_spec_files(log_files,avars.suffix[1])

            pgui('num_iq_files: ' + str(num_iq_files))
            pgui('num_log_files: ' + str(num_log_files))

            if(copy_extras == True):
                for j in range(len(avars.extra)):
                    avars.num_ex_files[j] = self.copy_spec_files(extra_files[j],avars.extra[j])
                pgui('num_ex_files = '+ str(avars.num_ex_files))     

            pgui('wrote %i sas files to %s\n' %
                        (num_iq_files, avars.sas_output_path))

            count += 1

        log.debug(vars(mvars))
        log.debug(vars(avars))

        return

    def epilogue(self):
        '''
        method to print out results and to move results
        to appropriate places.
        '''

        log = self.log
        pgui = self.run_utils.print_gui

        log.debug('in epilogue')

        self.run_utils.clean_up(log)

        pgui("%s \n" % ('=' * 60))
        pgui('EXTRACT UTILITIES IS DONE\n')
        pgui("%s \n" % ('=' * 60))
        time.sleep(1.0)

        return



