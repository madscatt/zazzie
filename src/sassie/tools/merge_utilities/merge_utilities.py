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
#	MERGE_UTILITIES
#
#	10/24/2012	--	initial coding			:	jc
#	04/19/2015	--	added dcd/pdb options   :	jc
#
#	 1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#								       *      **
'''
	MERGE_UTILITIES is a module that allows one to merge coordinates and
	scattering profiles from a series of runs.  Note that each coordinate
	scattering profile set is for the same molecule.

	INPUT: 
		Name of output directory		
		List of run_name directories and dcd filenames	
		List of scattering path names and sas profile type
		Extraction option

	OUTPUT:
	
		New single DCD file concatenated from multiple runs
		New single SAS profile directory with re-numbered filenames	
'''

import os
import sys
import locale
import numpy
import string
import time
import glob
import sassie.interface.input_filter as input_filter
import sasmol.sasmol as sasmol
import sassie.util.module_utilities as module_utilities
import sassie.util.sasconfig as sasconfig
import sassie.util.folder_management as folder_management


if sasconfig.__level__ == "DEBUG":
    DEBUG = True

app = 'merge_utilities'

class module_variables():

    def __init__(self, parent=None):
        self.app = app

class merge_utilities_input_variables():

    def __init__(self, parent=None):
        pass


class merge_utilities():

    def __init__(self, parent=None):
        pass

    def main(self, input_variables, txtOutput):

        self.mvars = module_variables()

        self.avars = merge_utilities_input_variables()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.run_utils.general_setup(self)

        self.initialization()

        self.merge_utilities()

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
        This method returns the variables to be used in the program.

       '''

        log = self.log
        mvars = self.mvars
        log.debug('in unpack_variables')

        mvars.runname = variables['runname'][0]
        mvars.pdb_file = variables['pdb_file'][0]                       #option (0==merge dcd/pdb/sas profiles, 1==merge sas only, 2==merge dcd/pdb only
        mvars.merge_option = variables['merge_option'][0]
        mvars.number_of_runs = variables['number_of_runs'][0]
        mvars.trajectory_names = variables['trajectory_names'][0]       # name of trajectory files to merge
        mvars.sas_type = variables['sas_type'][0]                       # sas profile type (0=sascalc,1=xtal2sas,2=cryson,3=crysol)
        mvars.sas_paths = variables['sas_paths'][0]                     # name of path with old sas profiles
        mvars.merge_type_option = variables['merge_type_option'][0]     # option (0==all, 1==weight files, 2==periodic)
        mvars.local_value = variables['local_value'][0]                 # None, list of weight files, or periodic value to skip
        mvars.output_filename = variables['output_filename'][0]         # name of output pdb or dcd file (if option 0 or 2)

        log.debug(vars(mvars))

        return 

    def copy_spec_files(self,sas_files, num_files, this_sas_output_path, suffix):

        log = self.log
        mvars = self.mvars
        avars = self.avars

        log.debug('in copy spec files')

        string_fill = 5
        number = 1
        log.debug('suffix: ' + suffix)

        for name in sas_files:
            if str(number) in avars.mask:
                numst = str(num_files)
                numst = numst.zfill(string_fill)
                runst = mvars.runname + '_' + numst + suffix[1:]
                this_file = os.path.join(this_sas_output_path, runst)
                cpst = 'cp  ' + name + ' ' + this_file
                log.debug('cpst: ' + cpst)
                os.system(cpst)
                num_files += 1

            number += 1

        return num_files


    def check_input_file_type(self,file_name):

        log = self.log
        log.debug ('in check input file type')
        
        input_type = None

        file_exist = os.path.isfile(file_name)

        if file_exist:
            mol = sasmol.SasMol(0)
            binary = input_filter.check_binary(file_name)

            if binary:
                mol.read_single_dcd_step(file_name, 0)
                input_type = 'dcd'
            else:
                mol.read_pdb(file_name, fastread=True)
                input_type = 'pdb'

        log.debug('input type: ' + input_type)
        
        return input_type


    def get_weight_file(self,local_value,**kwargs):

        log = self.log

        log.debug('in get weight files')
        log.debug('local_value: '+ str(local_value))
        
        print 'local_value = ', local_value

        weights = numpy.loadtxt(local_value)
        fweights = weights[:, 2]
        file_numbers = weights[:, 0]
        x2 = weights[:, 1]
        for i in xrange(len(x2)):
            if(fweights[i] == 1):
                log.debug('st = '+str(weights[i][0])+' : x2 = '+str(x2[i]))
                print 'st = ', int(weights[i][0]), ' : x2 = ', x2[i]

        try:
            coordinate_flag = kwargs['coordinate_flag']
#            print 'coordinate flag in get_weight_files: ',coordinate_flag
#            print 'merging coordinates: local_value = ', local_value
            coordinate_flag = True            
        except:
            coordinate_flag = False

        print 'coordinate flag after check: ', coordinate_flag
        
        if not coordinate_flag:
            frame_list = ((numpy.nonzero(fweights))[0]).tolist()
            mask = [str(int(file_numbers[i])) for i in frame_list]
            log.debug('mask (coord flag False): '+ str(mask))
        elif coordinate_flag:
            mask = ((numpy.nonzero(fweights))[0]).tolist()
            log.debug('mask (coord flag True): '+ str(mask))

#        print 'mask in get weight file: ', mask

        return mask


    def get_frequency(self, local_value,number_of_spec_files, **kwargs):

        log = self.log

        log.debug('in get frequency')
        log.debug('local_value: '+ str(local_value))
        log.debug('nspecfil: '+ str(number_of_spec_files))
        
        mask = []

#        print 'local value: ', local_value
#        print 'number of spec files: ', number_of_spec_files

        try:
            coordinate_flag = kwargs['coordinate_flag']
#            print 'coordinate flag in get_frequency: ',coordinate_flag
#            print 'merging coordinates: local_value = ', local_value
            coordinate_flag = True

        except:
            coordinate_flag = False

#        print 'coordinate flag after check: ', coordinate_flag

        if not coordinate_flag:
            for number in xrange(1, number_of_spec_files + 1, int(local_value)):
                mask.append(str(number))
        elif coordinate_flag:
            for number in xrange(0, number_of_spec_files, int(local_value)):
                mask.append(number)

        log.debug('mask in get frequency: '+ str(mask))
#        print 'mask in get frequency: ', mask

        return mask


    def get_sas_mask(self, local_value,**kwargs):

        log = self.log
        mvars = self.mvars

        log.debug('getting SAS mask')
        log.debug('merge type option (get mask): ' + str(mvars.merge_type_option))

#        print 'merge type option (get SAS mask): ', mvars.merge_type_option

        if mvars.merge_type_option == 1:
            mask = self.get_weight_file(local_value)
        elif mvars.merge_type_option == 2:
            number_of_spec_files = kwargs['number_of_spec_files']
            mask = self.get_frequency(local_value,number_of_spec_files)
        elif mvars.merge_type_option == 0:
            number_of_spec_files = kwargs['number_of_spec_files']
            local_value = '1'
            mask = self.get_frequency(local_value,number_of_spec_files)

        log.debug('mask in get SAS mask: '+ str(mask))
#        print 'mask in get SAS mask: ', mask

        return mask


    def initialization(self):
        '''
        method to prepare for merging of trajectories and/or SAS files
        '''

        log = self.log
        log.debug('in initialization')
        pgui = self.run_utils.print_gui

        mvars = self.mvars
        avars = self.avars

        log.debug(vars(mvars))
        log.debug(vars(avars))

        ''' directory and file preparation '''

        avars.output_path = mvars.runname + '/merge_utilities/'
        log.debug('output_path: '+avars.output_path)

#        version = 'version 0.2 : 04/19/15 : jc'
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


        avars.output_file = os.path.join(avars.output_path, mvars.output_filename)
#        print 'output_file: ',avars.output_file

#        print 'local_value: ', mvars.local_value        

        if mvars.merge_type_option == 1:
            avars.local_value = string.split(mvars.local_value, ',')
        else:
            avars.local_value = mvars.local_value

        print 'merge_type_option = ',mvars.merge_type_option, ; print type(mvars.merge_type_option)
        print 'local_value = ',avars.local_value, ; print type(avars.local_value)
            

        if (mvars.merge_option == 0 or mvars.merge_option == 2):

            if(mvars.output_filename[-4:] == ".dcd"):
                mergest = 'merging into dcd files \n'
                pgui(mergest)
                avars.output_type = 'dcd'
            elif(mvars.output_filename[-4:] == ".pdb"):
                mergest = 'merging into pdb files \n'
                pgui(mergest)
                avars.output_type = 'pdb'


        if (mvars.merge_option == 0 or mvars.merge_option ==1):

#           print 'sas_paths in merge_sas_files: ', mvars.sas_paths
            avars.sas_paths = [x.strip() for x in mvars.sas_paths.split(',')]
#           print 'sas_paths = ',avars.sas_paths
            avars.sas_input_paths = []
            for sas_path in avars.sas_paths:
                base_paths = []
#               print 'sas_path: ',sas_path
                new_paths = []
                for root, dirs, files in os.walk(sas_path, topdown=False):
                    for name in dirs:
#                       print(os.path.join(root, name))
                        new_paths.append(os.path.join(root, name))
                        base_paths.append(os.path.basename(os.path.normpath(os.path.join(root, name))))

#               print 'new_paths: ', new_paths
#               print 'base_paths: ', base_paths
                if base_paths == []:
                    base_paths.append(os.path.basename(os.path.normpath(sas_path)))
                    new_paths.append(sas_path)
#               print 'base_paths = ', base_paths
#               print 'new_paths = ', new_paths
                avars.sas_input_paths.append(new_paths)

                avars.sas_output_paths = []
        
                if(mvars.sas_type == 0):
                    for base_path in base_paths:
#                       print 'base_path: ', base_path       
                        avars.sas_output_paths.append(os.path.join(
                                avars.output_path, 'sascalc', base_path))
                    
#                   print 'sas_output_paths = ', avars.sas_output_paths
#                   print 'length sas_output_paths: ',len(avars.sas_output_paths)

                    avars.suffix = ['*.iq', '*.log']
                    avars.extra = []
                    avars.num_ex_files = [0,0,0,0]
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

            print 'sas_output_paths = ', avars.sas_output_paths
#           print 'length sas_output_paths: ',len(avars.sas_output_paths)
            print 'sas_input_paths = ', avars.sas_input_paths 
#           print 'length sas_input_paths: ',len(avars.sas_input_paths)                                            

        log.debug(vars(mvars))
        log.debug(vars(avars))


    def merge_utilities(self):
        '''
	    MERGE_UTILITIES is a module that allows one to merge coordinates and
	    scattering profiles from a series of runs.  Note that each coordinate
	    scattering profile set is for the same molecule.

        INPUT: variable descriptions

        runname:                string      project name
        pdb_file                string      input pdb file
        trajectory_names        string      input dcd files to be merged (number of files = number of runs)
        output_filename:        string      output dcd file 
        number_of_runs:         integer     number of dcd files and/or SAS runs to be merged                                           
        local_value:            string      value of merge_type_option (not used, list of weight file names, periodic value to skip) 
        merge_option:           integer     merge option (0==merge dcd/pdb/sas profiles, 1==merge sas only, 2==merge dcd/pdb only )                 
        merge_type_option:      integer     merge type option (0==all, 1==weight files, 2==periodic)               
        sas_type:               integer     integer depends on SAS file type (0==SasCalc, 1==Xtal2sas, 2==Cryson, 3==Crysol)    
        sas_paths:              string      paths to SAS files 

	    OUTPUT:
	
		New single DCD file concatenated from multiple runs 
		    and/or
		New single SAS profile directory with re-numbered filenames	

        '''

        log = self.log
        pgui = self.run_utils.print_gui
        log.debug('in merge_utilities')

        mvars = self.mvars
        avars = self.avars

        log.debug(vars(mvars))
        log.debug(vars(avars))

        ttxt = time.asctime(time.gmtime(time.time()))
        st = ''.join(['=' for x in xrange(60)])

        pgui("\n%s \n" % (st))
        pgui("DATA FROM RUN: %s \n\n" % (ttxt))

#       merge.log needs to go into sascalc directory so directories below sascalc can be selected in sassie-web
        if mvars.sas_type == 0:
            merge_log_path = avars.output_path+'/sascalc/'
            direxist = os.path.exists(merge_log_path)
            if(direxist == 0):
                try:
                    result = os.system('mkdir -p ' + merge_log_path)
                except:
                    message = 'can not create project directory: ' + merge_
                    log_path
                    message += '\nstopping here\n'
                    pgui(message)
                    sys.exit(1)
            avars.output_log_file = open(merge_log_path + 'merge.log', 'w')
        else:        
            avars.output_log_file = open(avars.output_path + 'merge.log', 'w')
        avars.output_log_file.write("DATA FROM RUN: %s \n\n" % (ttxt))


        if(mvars.merge_option == 0 or mvars.merge_option == 2):

            mergest = 'merging trajectory files\n'
            pgui(mergest)
            avars.trajectory_names = string.split(mvars.trajectory_names, ',')

            self.merge_trajectory_files()

            fraction_done = 0.5
            report_string = 'STATUS\t' + str(fraction_done)
            pgui(report_string)

        if(mvars.merge_option == 0 or mvars.merge_option == 1):

            mergest = '\n merging sas files\n'
            pgui(mergest)

            self.merge_sas_files()

        avars.output_log_file.close()

        fraction_done = 1.0
        report_string = 'STATUS\t' + str(fraction_done)
        pgui(report_string)

        return

    def merge_trajectory_files(self):

        log = self.log
        log.debug('in merge trajectory files')
        pgui = self.run_utils.print_gui

        mvars = self.mvars
        avars = self.avars

        log.debug(vars(mvars))
        log.debug(vars(avars))

        cpst = 'cp ' + mvars.pdb_file + ' ' + avars.output_path + os.path.sep
        print 'cpst: ',cpst
        os.system(cpst)

        m2 = sasmol.SasMol(1)
        m2.read_pdb(mvars.pdb_file, fastread=True)
        natoms = m2.natoms()

        if(avars.output_type == 'dcd'):
            dcdoutfile = m2.open_dcd_write(avars.output_file)

        m1 = sasmol.SasMol(0)
        m1.read_pdb(mvars.pdb_file, fastread=True)
        k = 0

        for i in xrange(len(avars.trajectory_names)):

#            print 'i, trajectory: ', i, avars.trajectory_names[i]
            avars.input_type = self.check_input_file_type(avars.trajectory_names[i])

#            print 'input type: ', avars.input_type

            if(avars.input_type == 'dcd'):
                dcdfile = m1.open_dcd_read(avars.trajectory_names[i])
                number_of_frames = dcdfile[2]
            if(avars.input_type == 'pdb'):
                pdbfile = m1.read_pdb(avars.trajectory_names[i])
                number_of_frames = m1.number_of_frames()


            avars.rangelist = [num for num in xrange(number_of_frames)]

            if(mvars.merge_type_option == 1):     # weight_files
                avars.rangelist = self.get_weight_file(avars.local_value[i], coordinate_flag=True)
                avars.step = 1
            elif(mvars.merge_type_option == 2):     # periodic
                avars.step = int(avars.local_value)
#                print 'local_value = ',avars.local_value, ; print type(avars.local_value)
            else:
                avars.step = 1
#            print 'step = ', avars.step

            pgui("reading %i frames from %s \n" %
                      (number_of_frames, avars.trajectory_names[i]))
            avars.output_log_file.write("reading %i frames from %s \n" %
                      (number_of_frames, avars.trajectory_names[i]))

            coor = numpy.zeros((1, natoms, 3), numpy.float32)


            for j in range(0, number_of_frames, avars.step):
                if j in avars.rangelist:
                    print '.',
                    sys.stdout.flush()
                    if(avars.input_type == 'dcd'):
                        m1.read_dcd_step(dcdfile, j)
                        coor[0, :, :] = m1.coor()[0]
                    elif(avars.input_type == 'pdb'):
                        coor[0, :, :] = m1.coor()[j]

                    m2.setCoor(coor)
                    if(avars.output_type == 'dcd'):
                        m2.write_dcd_step(dcdoutfile, 0, k + 1)
                    elif(avars.output_type == 'pdb'):
                        if(k == 0):
                            m2.write_pdb(avars.output_file, 0, 'w', model=k + 1)
                        else:
                            m2.write_pdb(avars.output_file, 0, 'a', model=k + 1)
                    k += 1

        pgui('\n')

        pgui('wrote %i frames to %s\n' % (k, './' + avars.output_file))
        avars.output_log_file.write('wrote %i frames to %s\n' % (k, './' + avars.output_file))

        if(avars.output_type == 'dcd'):
            m2.close_dcd_write(dcdoutfile)
        elif(avars.output_type == 'pdb'):
            with open(avars.output_file, "a") as myfile:
                myfile.write("END\n")

        log.debug(vars(mvars))
        log.debug(vars(avars))

        return


    def merge_sas_files(self):

        log = self.log
        log.debug('in merge sas files')
        pgui = self.run_utils.print_gui

        mvars = self.mvars
        avars = self.avars

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


        for i in xrange(len(avars.sas_output_paths)):
            copy_extras = False
            num_iq_files = 1
            num_log_files = 1
            for j in xrange(len(avars.sas_input_paths)):
                this_path = avars.sas_input_paths[j][i]
#               print 'i,j,this_path = ', i,j,this_path
                sys.stdout.flush()
                spec_files = []
                log_files = []
                extra_files = []
                for name in glob.glob(os.path.join(this_path, avars.suffix[0])):
                    spec_files.append(name)
                for name in glob.glob(os.path.join(this_path, avars.suffix[1])):
                    log_files.append(name)
                for ex in avars.extra:
                    this_extra = []
                    for name in glob.glob(os.path.join(this_path, ex)):
                        this_extra.append(name)
                    if(copy_extras == False and len(this_extra) > 0):
                        copy_extras = True
                        pgui('copying extra sas files')
                    this_extra.sort()
                    extra_files.append(this_extra)

                spec_files.sort()
                log_files.sort()

#                print 'spec_files: ', spec_files
#                print 'local_value: ',avars.local_value


                if(mvars.merge_type_option == 0):
                    avars.mask = self.get_sas_mask(avars.local_value,number_of_spec_files=len(spec_files))
                elif(mvars.merge_type_option == 1):
#                    print 'j,local_value[j]: ',j,avars.local_value[j]
                    avars.mask = self.get_sas_mask(avars.local_value[j])
                elif(mvars.merge_type_option == 2):
                    avars.mask = self.get_sas_mask(avars.local_value,number_of_spec_files=len(spec_files))

                num_iq_files = self.copy_spec_files(
                                    spec_files, num_iq_files, avars.sas_output_paths[i], avars.suffix[0])
                num_log_files = self.copy_spec_files(
                                    log_files, num_log_files, avars.sas_output_paths[i], avars.suffix[1])
                                
#               print 'i,num_iq_files: ', i, num_iq_files
                if(copy_extras == True):
                    for k in xrange(len(avars.extra)):
                        avars.num_ex_files[k] = self.copy_spec_files(extra_files[k], avars.num_ex_files[
                                                  k], avars.sas_output_paths[i], avars.extra[k])

            pgui('wrote %i sas files to %s\n' %
                            (num_iq_files - 1, './' + avars.sas_output_paths[i]))
#            print('wrote %i sas files to %s\n' %
#                            (num_iq_files - 1, './' + avars.sas_output_paths[i]))
            avars.output_log_file.write('wrote %i sas files to %s\n' %
                            (num_iq_files - 1, './' + avars.sas_output_paths[i]))

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

        st = ''.join(['=' for x in xrange(60)])
        pgui("\n%s \n" % (st))
        pgui('\nMERGE UTILITIES IS DONE')
        time.sleep(1.0)

        return



