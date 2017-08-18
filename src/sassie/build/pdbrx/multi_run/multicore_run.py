import os, glob
import multiprocessing 

import traceback
import time
import sys
sys.path.append('./')
import gui_mimic_pdbrx

class Process(multiprocessing.Process):
    def __init__(self, *args, **kwargs):
        multiprocessing.Process.__init__(self, *args, **kwargs)
        self._pconn, self._cconn = multiprocessing.Pipe()
        self._exception = None

    def run(self):
        try:
            multiprocessing.Process.run(self)
            self._cconn.send(None)
        except Exception as e:
            tb = traceback.format_exc()
            self._cconn.send((e, tb))
            # raise e  # You can still rise this exception if you need to

    @property
    def exception(self):
        if self._pconn.poll():
            self._exception = self._pconn.recv()
        return self._exception


def run_pdbrx(start, end, all_pdb_files, completed_pdb_files_directory, failed_pdb_files_directory):

	#time.sleep(10)
	fail = False
        test = False

	for job in xrange(end - start + 1):
		this_file = all_pdb_files[job + start]
		head, tail = os.path.split(this_file)	

                this_runname = 'run_'+tail[:-4]

                try:
                    run_gui = gui_mimic_pdbrx.gui_mimic_pdbrx(test, \
                                          runname=this_runname, \
                                          pdbfile=this_file, \
                                          use_defaults=True )
                    os.system('mv ' + this_runname + ' ' +
                              completed_pdb_files_directory)
                except:
                    os.system('mv ' + this_runname + ' ' +
                              failed_pdb_files_directory)
                              

		#print 'cp ' + this_file + ' ' + completed_pdb_files_directory + os.sep + tail + '_' + str(start)
		#os.system('cp ' + this_file + ' ' + completed_pdb_files_directory + os.sep + tail + '_' + str(start))

#		if job == 1:
#			fail = True	
		if fail:
			raise ValueError('Something went wrong...')

	return 


def main(all_pdb_files_directory, completed_pdb_files_directory, failed_pdb_files_directory, number_of_cores):
	
	all_pdb_files = glob.glob(all_pdb_files_directory + os.sep + '*.pdb')

	number_of_files = len(all_pdb_files)

	print 'number of pdb files found = ', number_of_files

	if(not os.path.exists(completed_pdb_files_directory)):
		os.system('mkdir ' + completed_pdb_files_directory)

	if(not os.path.exists(failed_pdb_files_directory)):
		os.system('mkdir ' + failed_pdb_files_directory)

	start = [] ; end = []

	remainder = number_of_files % number_of_cores
	divisor = number_of_files / number_of_cores
	print 'divisor = ', divisor
	print 'remainder = ', remainder

	'''
number of pdb files found =  7
divisor =  2
remainder =  1
start =  [0, 3, 6]
end =  [2, 5, 6]
	'''

        if divisor == 1:
            start.append(0)
            end.append(0)
        else:
	    for i in xrange(number_of_cores):
		    if i < number_of_cores - 1:
	                this_start = i * (divisor + 1)
	                start.append(this_start)
		        end.append(this_start + divisor)
	    else:
	        this_start = i * (divisor + 1)
	        start.append(this_start)
	        end.append(this_start + remainder - 1)

	print 'start = ', start
	print 'end = ', end

	all_jobs = []
        q = multiprocessing.Queue()

	for i in xrange(number_of_cores):
  
     		try: 
			this_process = multiprocessing.Process(target=run_pdbrx, args=(start[i], end[i], all_pdb_files, completed_pdb_files_directory, failed_pdb_files_directory))
			this_process.start()
		except Exception as error:
			print 'error = ', error
   
		all_jobs.append(this_process)

	return

if __name__ == '__main__':

	number_of_cores = 1

	all_pdb_files_directory = 'nr_pdb_files'
	completed_pdb_files_directory = 'completed_pdb_files'
	failed_pdb_files_directory = 'failed_pdb_files'

	main(all_pdb_files_directory, completed_pdb_files_directory, failed_pdb_files_directory, number_of_cores)
