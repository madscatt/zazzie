import os
import glob
import multiprocessing

import traceback
import time
import math
import sys
sys.path.append('./')
import gui_mimic_pdbrx

def find_and_move_json_log(this_file, this_runname, failed_pdb_files_directory):

    all_json_files = glob.glob('*_json')
    all_log_files = glob.glob('*_log')

    for this_json in all_json_files:
        result = open(this_json).readlines()
        print 'this_json = ', this_json
        print 'result = ', result
        for line in result:
            if this_file in line:
                mvst = 'mv -f ' + this_json + ' ' + failed_pdb_files_directory + \
                os.sep + this_runname + os.sep + 'pdbrx' + os.sep
                print mvst
                print mvst
                print mvst
                os.system(mvst)
                break 

    for this_log in all_log_files:
        result = open(this_log).readlines()
        for line in result:
            if this_file in line:
                mvst = 'mv -f ' + this_log + ' ' + failed_pdb_files_directory + \
                os.sep + this_runname + os.sep + 'pdbrx' + os.sep
                print mvst
                print mvst
                print mvst
                os.system(mvst)
                break 

    return

def run_pdbrx(start, end, all_pdb_files, completed_pdb_files_directory, failed_pdb_files_directory):

        # time.sleep(10)
    fail = False
    test = False

    for job in xrange(end - start + 1):
        this_file = all_pdb_files[job + start]
        head, tail = os.path.split(this_file)

        this_runname = 'run_' + tail[:-4]

        try:
            run_gui = gui_mimic_pdbrx.gui_mimic_pdbrx(test,
                                                      runname=this_runname,
                                                      pdbfile=this_file,
                                                      use_defaults=True)
            os.system('mv -f ' + this_runname + ' ' +
                      completed_pdb_files_directory)
            os.system('mv -f ' + this_file + ' ' +
                      completed_pdb_files_directory)
        except:
            os.system('mv -f ' + this_runname + ' ' +
                      failed_pdb_files_directory)
           
            find_and_move_json_log(this_file, this_runname, failed_pdb_files_directory) 


        time.sleep(1)

        # print 'cp ' + this_file + ' ' + completed_pdb_files_directory + os.sep + tail + '_' + str(start)
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

    start, end, number = get_ranges(number_of_cores, number_of_files)

    print '='*50
    print 'start = ', start
    print 'end = ', end
    print 'number = ', number
    print '='*50

    #sys.exit()

    all_jobs = []
    q = multiprocessing.Queue()

    for i in xrange(number_of_cores):

        try:
            this_process = multiprocessing.Process(target=run_pdbrx, args=(
                start[i], end[i], all_pdb_files, completed_pdb_files_directory, failed_pdb_files_directory))
            this_process.start()
            time.sleep(1)

        except Exception as error:
            print 'error = ', error

        all_jobs.append(this_process)

    return


def get_ranges(ncpus, nfiles):

    if nfiles < ncpus:
        nbatches = 0
        nfiles_per_batch = nfiles
    else:
        nbatches = int(math.floor(float(nfiles) / float(ncpus)))
        nfiles_per_batch = ncpus

    #nfiles_per_batch = nfiles / nbatches

    if nbatches != 1 :
        remainder = nfiles % ncpus
    else:
        remainder = 0

    remainder = nfiles % ncpus

    print 'ncpus = ', ncpus
    print 'nfiles = ', nfiles
    print 'nbatches = ', nbatches
    print 'nfiles_per_batch = ', nfiles_per_batch
    print 'remainder = ', remainder

    first = [] ; last = [] ; number = []

    for i in xrange(nbatches):
        this_first = i * nfiles_per_batch
        first.append(this_first)
        this_last = this_first + nfiles_per_batch - 1
        last.append(this_last)
        number.append(this_last - this_first + 1)
        print 'i = ', i, ' : first = ', this_first, ' : last = ', this_last, ' : number = ', number[-1]

    if nbatches > 1: 
        final_first = this_last + 1
        final_last = final_first + remainder - 1
        final_number = final_last - final_first + 1
        print 'final : first = ', final_first, ' : last = ', final_last, ' : number = ', final_number

    else:
        final_first = 0
        final_last = remainder - 1
        final_number = remainder
        print 'final : first = ', final_first, ' : last = ', final_last, ' : number = ', final_number

    if final_number > 0:
        first.append(final_first)
        last.append(final_last)
        number.append(final_number)

    return first, last, number


if __name__ == '__main__':

    number_of_cores = 64

    all_pdb_files_directory = 'nr_pdb_files'
    completed_pdb_files_directory = 'completed_pdb_files'
    failed_pdb_files_directory = 'failed_pdb_files'

    main(all_pdb_files_directory, completed_pdb_files_directory,
         failed_pdb_files_directory, number_of_cores)

