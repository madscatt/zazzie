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
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
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

import os,sys,locale,numpy,string,time,glob
import sassie.interface.input_filter as input_filter
try:
    import sassie.sasmol.sasmol as sasmol
except:
    import sasmol.sasmol as sasmol
try:
    import sassie.sasconfig as sasconfig
except:
    import sassie.util.sasconfig as sasconfig

sys.path.append('./')

def print_failure(message,txtOutput):

        txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
        txtOutput.put(">>>> RUN FAILURE <<<<\n")
        txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
        txtOutput.put(message)

        return

def unpack_variables(variables):
	'''
	This method returns the variables to be used in the program.

	'''

	runname		= variables['runname'][0]
	pdb_file	= variables['pdb_file'][0]
	merge_option	= variables['merge_option'][0]	# option (0==merge dcd/pdb/sas profiles, 1==merge sas only, 2==merge dcd/pdb only 
	number_of_runs 	= variables['number_of_runs'][0] 	
	trajectory_names 	= variables['trajectory_names'][0]	# name of trajectory files to merge
	sas_type 	= variables['sas_type'][0] 	# sas profile type (1=xtal2sas,2=cryson,3=crysol)
	sas_paths	= variables['sas_paths'][0]	# name of path with old sas profiles
	merge_type_option	= variables['merge_type_option'][0]	# option (0==all, 1==weight files, 2==periodic)
	local_value	= variables['local_value'][0]	# None, list of weight files, or periodic value to skip

	output_filename	= variables['output_filename'][0]	# name of output pdb or dcd file (if option 0 or 2)


	return runname,pdb_file,merge_option,number_of_runs,trajectory_names,sas_type,sas_paths,output_filename,merge_type_option,local_value


def check_input_file_type(file_name):

    input_type = None

    file_exist = os.path.isfile(file_name)

    if file_exist:
        mol = sasmol.SasMol(0)
        binary = input_filter.check_binary(file_name)

        if binary:
            mol.read_single_dcd_step(file_name,0)
            input_type = 'dcd'
        else:
            mol.read_pdb(file_name,fastread=True)     
            input_type = 'pdb'

    return input_type 


def get_weight_file(local_value, **kwargs):
    print 'local_value = ',local_value

    weights = numpy.loadtxt(local_value)
    fweights = weights[:,2]
    file_numbers = weights[:,0]
    rg = weights[:,1]
    for i in xrange(len(rg)):
        if(fweights[i]==1):
            print 'st = ',int(weights[i][0]),' : rg = ',rg[i]
      
    try:
        coordinate_flag = kwargs['coordinate_flag']
    except:
        coordinate_flag = False

    if not coordinate_flag:
        frame_list = ((numpy.nonzero(fweights))[0]).tolist()
        #mask = [str(i) for i in frame_list]
        mask = [str(int(file_numbers[i])) for i in frame_list]
    elif coordinate_flag:
        mask = ((numpy.nonzero(fweights))[0]).tolist()
   
    return mask

def get_frequency(local_value,number_of_spec_files, **kwargs):
  
    mask = []
     
    try:
        coordinate_flag = kwargs['coordinate_flag']
    except:
        coordinate_flag = False

    if not coordinate_flag:
        for number in xrange(1,number_of_spec_files+1,int(local_value)):
            mask.append(str(number))
    elif coordinate_flag:
        for number in xrange(0,number_of_spec_files,int(local_value)):
            mask.append(number)

    return mask

def get_sas_mask(extract_option, value, **kwargs):

    if extract_option == 'weight_file':
        mask = get_weight_file(value)
    elif extract_option == 'sampling_frequency':
        number_of_spec_files = kwargs['number_of_spec_files']
        mask = get_frequency(value,number_of_spec_files)
    elif extract_option == 'all':
        number_of_spec_files = kwargs['number_of_spec_files']
        value = '1'
        mask = get_frequency(value,number_of_spec_files)
     
    return mask 

def merge_trajectory_files(runname,pdb_file,trajectory_files,output_path,output_log_file,output_type,output_filename,merge_type_option,local_value,txtOutput):

    direxist=os.path.exists(output_path)
    if(direxist==0):
        try:
            result=os.system('mkdir -p '+output_path)
        except:
            message='can not create project directory: '+output_path
            message+='\nstopping here\n'
            print_failure(message,txtOutput)

    output_file = os.path.join(output_path,output_filename)

    cpst = 'cp '+pdb_file+' '+output_path+os.path.sep
    os.system(cpst)

    #txtOutput.put("writing frames to %s \n" % (output_file))
    #output_log_file.write("writing frames to %s \n" % (output_file))
	
    m2 = sasmol.SasMol(1)	
    m2.read_pdb(pdb_file,fastread=True)
    natoms = m2.natoms()	

    if(output_type == 'dcd'):
        dcdoutfile = m2.open_dcd_write(output_file)

    m1 = sasmol.SasMol(0)
    m1.read_pdb(pdb_file,fastread=True)
	
    k=0
    for i in xrange(len(trajectory_files)):

        input_type = check_input_file_type(trajectory_files[i])

        if(input_type == 'dcd'):
            dcdfile = m1.open_dcd_read(trajectory_files[i])
            number_of_frames = dcdfile[2]
        if(input_type == 'pdb'):
            pdbfile = m1.read_pdb(trajectory_files[i])
            number_of_frames = m1.number_of_frames()
    
        #print 'merge_type_option = ',merge_type_option, ; print type(merge_type_option)
        #print 'local_value = ',local_value, ; print type(local_value)

        rangelist = [num for num in xrange(number_of_frames)]

        if(merge_type_option == 1):     # weight_files
            rangelist = get_weight_file(local_value[i],coordinate_flag=True)
            step = 1
        elif(merge_type_option == 2):     # periodic
            step = int(local_value)
        #    print 'local_value = ',local_value, ; print type(local_value)
        else:
            step = 1

        print 'step = ',step

        txtOutput.put("reading %i frames from %s \n" % (number_of_frames,trajectory_files[i]))
        output_log_file.write("reading %i frames from %s \n" % (number_of_frames,trajectory_files[i]))
        coor = numpy.zeros((1,natoms,3),numpy.float32)
        for j in range(0,number_of_frames,step):
            if j in rangelist:
                print '.', ; sys.stdout.flush()
                if(input_type == 'dcd'):
                    m1.read_dcd_step(dcdfile,j)
                    coor[0,:,:] = m1.coor()[0]
                elif(input_type == 'pdb'):
                    coor[0,:,:] = m1.coor()[j]

                m2.setCoor(coor)
                if(output_type == 'dcd'):
                    m2.write_dcd_step(dcdoutfile,0,k+1)
                elif(output_type == 'pdb'):
                    if(k == 0):
                        m2.write_pdb(output_file,0,'w',model=k+1)
                    else:
                        m2.write_pdb(output_file,0,'a',model=k+1)
                k += 1
   
    print '\n'

    txtOutput.put('wrote %i frames to %s\n' % (k,'./'+output_file))
    output_log_file.write('wrote %i frames to %s\n' % (k,'./'+output_file))

    if(output_type == 'dcd'):
        m2.close_dcd_write(dcdoutfile)
    elif(output_type == 'pdb'):
        with open(output_file, "a") as myfile:
            myfile.write("END\n")

    return

def copy_spec_files(sas_files,num_files,runname,sas_output_path,suffix,mask):
	
    string_fill = 5
   
    number = 1
     
    for name in sas_files:
        if str(number) in mask:
            numst = str(num_files) ; numst=numst.zfill(string_fill)
            runst = runname+'_'+numst+suffix[1:]	
            this_file = os.path.join(sas_output_path,runst)
            cpst = 'cp  '+name+' '+this_file
            os.system(cpst)
            num_files+=1
        number += 1

    return num_files

def merge_sas_files(runname,all_sas_paths,sas_type,output_path,output_log_file,merge_type_option,local_value,txtOutput):

    if(sas_type==1):
        sas_output_path = output_path+'xtal2sas'
        suffix=['*.iq','*.log']
        extra = ['*.inf','*.crd','*.ans','*.pr']
        num_ex_files=[1,1,1,1]
    elif(sas_type==2):
        sas_output_path = output_path+'cryson'
        suffix=['*.int','*.log']
        extra = ['*.sav','*.flm','*.alm','*.ans']
        num_ex_files=[1,1,1,1]
    elif(sas_type==3):
        sas_output_path = output_path+'crysol'
        suffix=['*.int','*.log']
        extra = ['*.sav','*.flm','*.alm','*.ans']
        num_ex_files=[1,1,1,1]

    direxist=os.path.exists(sas_output_path)
    if(direxist==0):
        try:
            result=os.system('mkdir -p '+sas_output_path)
        except:
            message='can not create project directory: '+sas_output_path
            message+='\nstopping here\n'
            print_failure(message,txtOutput)

    copy_extras=False	

    num_iq_files = 1
    num_log_files = 1
	
    for i in xrange(len(all_sas_paths)):
        this_path = all_sas_paths[i]
        print 'this_path = ',this_path ; sys.stdout.flush()
        spec_files = [] ; log_files = []
        extra_files = []
        for name in glob.glob(os.path.join(this_path,suffix[0])):
            spec_files.append(name)		
        for name in glob.glob(os.path.join(this_path,suffix[1])):
            log_files.append(name)		
        for ex in extra:
            this_extra = []
            for name in glob.glob(os.path.join(this_path,ex)):
                this_extra.append(name)
            if(copy_extras==False and len(this_extra)>0):
                copy_extras=True
                print 'copying extra sas files'
            this_extra.sort()		
            extra_files.append(this_extra)	
	
        spec_files.sort()
        log_files.sort()

        if(merge_type_option == 0):
            extract_option = 'all'
            mask = get_sas_mask(extract_option, local_value,number_of_spec_files=len(spec_files))
        elif(merge_type_option == 1):
            extract_option = 'weight_file'
            mask = get_sas_mask(extract_option, local_value[i])
        elif(merge_type_option == 2):
            extract_option = 'sampling_frequency'
            mask = get_sas_mask(extract_option, local_value,number_of_spec_files=len(spec_files))
            	
        num_iq_files=copy_spec_files(spec_files,num_iq_files,runname,sas_output_path,suffix[0],mask)
        num_log_files=copy_spec_files(log_files,num_log_files,runname,sas_output_path,suffix[1],mask)

        if(copy_extras == True):		
            for j in xrange(len(extra)):
                num_ex_files[j]=copy_spec_files(extra_files[j],num_ex_files[j],runname,sas_output_path,extra[j],mask)
			
    txtOutput.put('wrote %i sas files to %s\n' % (num_iq_files-1,'./'+sas_output_path))
    output_log_file.write('wrote %i sas files to %s\n' % (num_iq_files-1,'./'+sas_output_path))

    return

def merge_runs(variables,txtOutput):

    runname,pdb_file,merge_option,number_of_runs,trajectory_names,sas_type,sas_paths,output_filename,merge_type_option,local_value = unpack_variables(variables)

    #ttxt=time.ctime()
    ttxt=time.asctime( time.gmtime( time.time() ) ) 
    st=''.join(['=' for x in xrange(60)])

    txtOutput.put("\n%s \n" %(st))
    txtOutput.put("DATA FROM RUN: %s \n\n" %(ttxt))

    output_path=runname+'/merge_utilities/' 
    print 'output_path = ',output_path

    version='version 0.2 : 04/19/15 : jc'
    direxist=os.path.exists(output_path)
    if(direxist==0):
        try:
            result=os.system('mkdir -p '+ output_path)
        except:
            message='can not create project directory: '+output_path
            message+='\nstopping here\n'
            print_failure(message,txtOutput)

    output_log_file = open(output_path+'merge.log','w')
    output_log_file.write("DATA FROM RUN: %s \n\n" %(ttxt))

    if merge_type_option == 1:
        local_value = string.split(local_value, ',')

    if(merge_option == 0 or merge_option == 2):

        if(output_filename[-4:] == ".dcd"):
            mergest =  '> merging into dcd files \n'
            print mergest
            output_type = 'dcd'
        elif(output_filename[-4:] == ".pdb"):
            mergest = '> merging into pdb files \n'
            print mergest
            output_type = 'pdb'

        mergest = '> merging trajectory files\n'
        txtOutput.put(mergest)
        trajectory_names = string.split(trajectory_names,',')

        merge_trajectory_files(runname,pdb_file,trajectory_names,output_path,output_log_file,output_type,output_filename,merge_type_option,local_value,txtOutput)

        fraction_done = 0.5
        report_string='STATUS\t'+str(fraction_done)
        txtOutput.put(report_string)
	
    if(merge_option == 0 or merge_option == 1):
             
        mergest = '\n> merging sas files\n'
        print mergest
        txtOutput.put(mergest)
        print 'BEFORE MERGE SAS'
        print 'sas_paths = ',sas_paths
        print 'type(sas_paths) = ',type(sas_paths)
        sas_paths = string.split(sas_paths,',')

        merge_sas_files(runname,sas_paths,sas_type,output_path,output_log_file,merge_type_option,local_value,txtOutput)
	
    output_log_file.close()

#
#       write global run name, pdb, and dcd filenames to .last_sas
#       
    final_dcd_file_name = runname+'.dcd'
    if(sas_type==1): 
        xtalpath = runname+'/merge/xtal2sas'
    elif(sas_type==2):
        xtalpath = runname+'/merge/cryson'
    elif(sas_type==3):
        xtalpath = runname+'/merge/crysol'

    fileexist=os.path.exists('.last_sas')
    if(fileexist==1):
        os.system('mv -f .last_sas .last_sas_bu')
    lastsasfile=open('./.last_sas','w')
    lastsasfile.write('runname\t'+runname+'\n')
    lastsasfile.write('pdb_name\t'+pdb_file+'\n')
    lastsasfile.write('dcd_name\t'+final_dcd_file_name+'\n')

    lastsasfile.write('sas_name\t'+xtalpath+'\n')
    lastsasfile.close()

    fraction_done = 1.0
    report_string='STATUS\t'+str(fraction_done)
    txtOutput.put(report_string)

    txtOutput.put("\n%s \n" %(st))
    print '\nMERGE UTILITIES IS DONE'
    time.sleep(2.5)

    return

def report_error(error):
    print 'error = ',error
    print 'QUITTING NOW'
    sys.exit()

if __name__ == '__main__':

    ###
    ###  BEGIN USER EDIT    ###
    ###

    runname = 'run_0'
    pdb_file = 'min3.pdb' 
    merge_option = '0'   # 0 -> pdb/dcd:sas ; 1 -> sas ; 2 -> pdb/dcd
    merge_option = '2'   # 0 -> pdb/dcd:sas ; 1 -> sas ; 2 -> pdb/dcd
    number_of_runs = '2' 
    trajectory_names = 'run_m1/generate/run_m1.dcd,run_m2/generate/run_m2.dcd'
    trajectory_names = 'min3.pdb,min3.pdb'
    trajectory_names = 'min3.pdb,run_m2/generate/run_m2.dcd'
    sas_type = '1'      # 1 -> xtal2sas, 2 -> cryson, 3 -> cyrsol
    sas_paths = 'run_m1/xtal2sas,run_m2/xtal2sas'
    merge_type_option = '0'   # 0 -> all, 1 -> weight files, 2 -> periodic
    #merge_type_option = '1'   # 0 -> all, 1 -> weight files, 2 -> periodic
    #merge_type_option = '2'   # 0 -> all, 1 -> weight files, 2 -> periodic
    local_value = 'None'
    #local_value = 'weights_file_m1.txt,weights_file_m2.txt'
    #local_value = '10'
    
    output_filename = 'merged_run_0.pdb' 
     
    ###
    ###  END USER EDIT      ###
    ###
    
    svariables = {} 
    svariables['runname'] = (runname,'string')
    svariables['pdb_file'] = (pdb_file,'string')
    svariables['merge_option'] = (merge_option,'int')
    svariables['number_of_runs'] = (number_of_runs,'int')
    svariables['trajectory_names'] = (trajectory_names,'string')
    svariables['sas_type'] = (sas_type,'int')
    svariables['sas_paths'] = (sas_paths,'string')
    svariables['merge_type_option'] = (merge_type_option,'int')
    svariables['local_value'] = (local_value,'string')
    
    svariables['output_filename'] = (output_filename,'string')
    
    import sassie.interface.input_filter as input_filter
    import sassie.interface.merge_utilities_filter as merge_utilities_filter   
    import multiprocessing

    print 'running as a main process'
    
    error,variables=input_filter.type_check_and_convert(svariables)
    error=merge_utilities_filter.check_merge_utilities(variables)

    if(len(error) != 0): report_error(error)

    txtOutput=multiprocessing.JoinableQueue()

    merge_runs(variables,txtOutput)

