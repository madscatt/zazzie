import os
import glob
import numpy
import locale
import string

DEBUG = True

if DEBUG:
    pass
#    debug_number_of_files = 100
#    debug_number_of_files = 500
#    debug_number_of_files = 10000

#TODO: need to make sure API on file reading is followed exactly
#TODO: need to have a policy for lines in data file that should be okay to skip
#TODO: need to look at UTF-8 encoding 
#TODO: need to handle new-line / carriage return 
#TODO: need example MSDOS/WINDOWS file example as a test case

#TODO: all the above should be captured in an API ; perhaps sassie-wide file reader with options

#TODO: these same file reading & testing methods should be used by the filter for this module

def validate_data(other_self):

    mvars = other_self.mvars
    evars = other_self.evars
    log = other_self.log
    pgui = other_self.run_utils.print_gui 

    if DEBUG:
        pgui('in validate_data')

        log.debug('q_data = ' + numpy.array2string(evars.q_data))
        log.debug('goal_q_data = ' + numpy.array2string(evars.goal_q_data))
        log.debug('goal_iq_data = ' + numpy.array2string(evars.goal_iq_data))
        log.debug('goal_iq_error_data = ' + numpy.array2string(evars.goal_iq_error_data))

    #TODO: use something like this to create an error: q_data[0] = 3.3333
    #TODO: all of this should be handled in filter

    #check that q-arrays of theoretical and experimental data are identical

    if evars.q_data.all() != evars.goal_q_data.all():
        pgui("q-values in theoretical files do not agree with those in experimental data file")

    #check lengths of theoretical arrays: number of frames
   
    #TODO: not tested for failure
    if evars.iq_data.shape != evars.iq_error_data.shape:
        pgui("number of theoretical iq values does not match theoretical iq error values")

    # check that the number of data points match
    
    #TODO: not tested for failure
    if evars.iq_data.shape[1] != evars.goal_q_data.shape[0]:
        pgui("number of theoretical iq values does not match the number of experimental q-values")
        pgui("evars.iq_data.shape[1] = "+str(evars.iq_data.shape[1]))
        pgui("evars.goal_q_data.shape[0] = "+str(evars.goal_q_data.shape[0]))

    #TODO: not tested for failure
    if evars.iq_error_data.shape[1] != evars.goal_q_data.shape[0]:
        pgui("number of theoretical iq error values does not match the number of experimental q-values")
        pgui("evars.iq_error_data.shape[1] = "+str(evars.iq_error_data.shape[1]))
        pgui("evars.goal_q_data.shape[0] = "+str(evars.goal_q_data.shape[0]))

    #TODO: not tested for failure
    if evars.goal_iq_data.shape[0] != evars.goal_q_data.shape[0]:
        pgui("number of experimental iq values does not match the number of experimental q-values")
        pgui("evars.goal_iq_data.shape[0] = "+str(evars.goal_iq_data.shape[0]))
        pgui("evars.goal_q_data.shape[0] = "+str(evars.goal_q_data.shape[0]))

    #TODO: not tested for failure
    if evars.goal_iq_error_data.shape[0] != evars.goal_q_data.shape[0]:
        pgui("number of experimental iq error values does not match the number of experimental q-values")
        pgui("evars.goal_iq_error_data.shape[0] = "+str(evars.goal_iq_error_data.shape[0]))
        pgui("evars.goal_q_data.shape[0] = "+str(evars.goal_q_data.shape[0]))


    return

def read_goal_iq_data(other_self):

    mvars = other_self.mvars
    evars = other_self.evars
    log = other_self.log
    pgui = other_self.run_utils.print_gui 

    this_q = []
    this_iq = []
    this_iq_error = []
    with open(mvars.goal_iq_data_file) as this_file:
        for line in this_file:  
            my_line = string.split(line)
            this_q.append(locale.atof(my_line[0]))
            this_iq.append(locale.atof(my_line[1]))
            this_iq_error.append(locale.atof(my_line[2]))

    this_q = numpy.array(this_q, numpy.float)
    this_iq = numpy.array(this_iq, numpy.float)
    this_iq_error = numpy.array(this_iq_error, numpy.float)

    evars.goal_q_data = this_q
    evars.goal_iq_data = this_iq
    evars.goal_iq_error_data = this_iq_error

    return 

def read_iq_data(other_self):
    mvars = other_self.mvars
    evars = other_self.evars
    log = other_self.log
    pgui = other_self.run_utils.print_gui 

    if DEBUG:
        pgui("DEBUG MODE")
        pgui("DEBUG MODE")
        pgui("DEBUG MODE")

    first = True

    evars.number_iq_files = 0

    #for name in glob.glob(os.path.join(mvars.iq_data_path, '*.iq')):    
    for name in glob.glob(os.path.join(mvars.iq_data_path, '*.ciq')):    
        evars.number_iq_files += 1
    pgui('evars.number_iq_files = ' + str(evars.number_iq_files))

    if mvars.number_of_files_to_use < evars.number_iq_files:
        log.debug("fewer files found than requested number of files by user")
        evars.number_iq_files = mvars.number_of_files_to_use
    if mvars.number_of_files_to_use > evars.number_iq_files:
        log.error("more files requested by user than found")

    count = 0

    #for name in glob.glob(os.path.join(mvars.iq_data_path, '*.ciq')):    
    for name in glob.glob(os.path.join(mvars.iq_data_path, '*.ciq')):    

        this_q = []
        this_iq = []
        this_iq_error = []
        with open(name) as this_file:
            for line in this_file:  
                my_line = string.split(line)
                this_q.append(locale.atof(my_line[0]))
                this_iq.append(locale.atof(my_line[1]))
                this_iq_error.append(locale.atof(my_line[2]))

        this_q = numpy.array(this_q, numpy.float)
        this_iq = numpy.array(this_iq, numpy.float)
        this_iq_error = numpy.array(this_iq_error, numpy.float)
        
        if first:

            if DEBUG:
                debug_number_of_files = mvars.number_of_files_to_use
                q_data = numpy.zeros(len(this_q), numpy.float)
                iq_data = numpy.zeros((debug_number_of_files,len(this_q)), numpy.float)
                iq_data = numpy.zeros((debug_number_of_files,len(this_q)), numpy.float)
                iq_error_data = numpy.zeros((debug_number_of_files,len(this_q)), numpy.float)
            else:
                q_data = numpy.zeros(len(this_q), numpy.float)
                iq_data = numpy.zeros((mvars.number_of_files_to_use, len(this_q)), numpy.float)
                iq_error_data = numpy.zeros((mvars.number_of_files_to_use, len(this_q)), numpy.float)

            q_data = this_q
            first = False        
        
        iq_data[count,:] = this_iq
        iq_error_data[count,:] = this_iq_error

        if count > 0:
            #test if all q-values are identical
            #if(count == 3):
            #    this_q[0] = 3.0
            if this_q.all() != last_q.all():
                pgui('q-values do not agree for file: '+name+'\nSTOPPING PROGRAM\n')
                log.error('q-values do not agree for file: '+name+'\nSTOPPING PROGRAM\n')
                pgui('this_q = ' + str(this_q))
                log.error('this_q = ' + str(this_q))
                pgui('last_q = ' + str(last_q))
                log.error('last_q = ' + str(last_q))
                return False

        last_q = this_q
             
        count += 1
        if(count > (debug_number_of_files - 1) and DEBUG):
            log.debug('q_data TEST final q = ' + numpy.array2string(q_data[-1]))
            log.debug('iq_data TEST final q = ' + numpy.array2string(iq_data[count-1,-1]))
            log.debug('iq_error_data TEST final q = ' + numpy.array2string(iq_error_data[count-1,-1]))
            evars.q_data = q_data
            evars.iq_data = iq_data
            evars.iq_error_data = iq_error_data
            evars.number_iq_files = debug_number_of_files

            return 

        evars.q_data = q_data
        evars.iq_data = iq_data
        evars.iq_error_data = iq_error_data

    return 
