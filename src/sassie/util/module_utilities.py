import sys
import os
import logging
import getpass
import platform
import time
#import pkg_resources
import json


import sassie.util.sasconfig as sasconfig

#import sasconfig 

if sasconfig.__level__ == "DEBUG": DEBUG = True

class run_utils():

    def __init__(self,app,txtOutput):
        self.__application__ = app
        self.txtOutput = txtOutput

    def write_json_file(self):

        self.logger.debug('in write_json_file')

        self.logger.info('writing json files with input variables')

        with open(self.parameter_file,'w') as outfile:
            json.dump([self.v],outfile)
            self.logger.info('writing parameters to '+self.parameter_file)

        if(os.path.isfile('.last_sas.json')):
            os.system('mv .last_sas.json .last_sas.json_backup')
            self.logger.info('backing up existing .last_sas.json file to .last_sas.json_backup')

        if(os.path.isfile('.last_sas_'+self.__application__+'.json')):
            os.system('mv .last_sas_'+self.__application__+'.json .last_sas_'+self.__application__+'.json_backup')
            self.logger.info('backing up existing .last_sas_'+self.__application__+'.json file to .last_sas_'+self.__application__+'.json_backup')

        with open('.last_sas.json','w') as outfile:
            json.dump([self.v],outfile)
            self.logger.info('writing parameters to .last_sas.json file')

        with open('.last_sas_'+self.__application__+'.json','w') as outfile:
            json.dump([self.v],outfile)
            self.logger.info('writing parameters to .last_sas_'+self.__application__+'.json file')

        self.logger.info('input variables: '+json.dumps([self.v] ))

        return

    def general_setup(self,other_self):

        ''' method to write json file of input variables and setup logging '''

        ''' grab all of the variables in the variable class instance '''
        log = other_self.log
        module_variables = other_self.module_variables

        log.debug('in general_setup')

        input_vars = [attr for attr in dir(module_variables) if not callable(getattr(module_variables,attr)) and not attr.startswith("__")]

        ''' put the variables into a dictionary to send to json '''

        input_variables = {}

        for var in input_vars:
            input_variables[var] = getattr(module_variables, var)

        self.v = input_variables

        self.write_json_file()
        run_name = self.v['run_name']
        self.run_path = os.path.join(run_name,self.__application__)

        self.logger.info('setting run_path directory to : '+self.run_path)

        try:
            if(not os.path.exists(run_name)):
                os.makedirs(self.run_path)
            elif(not os.path.exists(self.run_path)):
                os.mkdir(self.run_path)
            self.logger.info('creating directory : '+self.run_path)
        except:
            self.logger.critical('FAILED to create '+self.run_path+' directory')

        other_self.run_path = self.run_path
        other_self.parameter_file = self.parameter_file

        return

    def preamble(self):

        user = getpass.getuser()
        current_path = os.getcwd()
        self.logger.debug('in preamble')
        #self.logger.info('sassie version(sassie) : '+sassie.__version__)
        #self.logger.info('sassie version : '+pkg_resources.get_distribution("sassie_2").version)
        self.logger.info('module_application : '+self.__application__)
        self.logger.info('executed by user : '+user+' on '+time.ctime())
        self.logger.info('current directory : '+current_path)
        self.logger.info('hostname : '+platform.node())

        return

    def setup_logging(self,other_self):

        self.logger = logging.getLogger(self.__application__)
        if (sasconfig.__level__ == 'DEBUG'):
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.INFO)

        time_string = time.strftime("%Y%m%d-%H%M%S")

        self.log_file = self.__application__+'_'+time_string+'.sassie_log'
        self.parameter_file = self.log_file[:-3]+'json'
        outfile = open(self.log_file,"w")
        st = 'Date\t\tTime\t\tFile\t\tMethod\t\tLine\tLevel\tMessage\n'
        outfile.write(st) #; print(st)
        outfile.close()

        file_handler = logging.FileHandler(self.log_file)
        formatter = logging.Formatter('%(asctime)s - %(filename)s - %(funcName)s - %(lineno)d - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)

        output_handler = logging.StreamHandler(sys.stdout)
        output_handler.setLevel(logging.WARNING)
        output_handler.setFormatter(formatter)

        self.logger.addHandler(file_handler)
        self.logger.addHandler(output_handler)

        self.logger.debug('in setup_logging')
        self.logger.info("Program started")

        self.preamble()

        other_self.log = self.logger
        other_self.log_file = self.log_file

        return

    def print_gui(self,message):
        '''
        method to output message to logger and txtOutput (gui)
        '''

        self.logger.info(message)
        self.txtOutput.put(message)
        if DEBUG: print(message)

    def clean_up(self,log):
        '''
        method to move files to run_path directory and finish tasks
        '''

        log.debug('in clean_up')

        os.system('mv '+self.log_file+' '+self.run_path)
        os.system('mv '+self.parameter_file+' '+self.run_path)

    def capture_exception(self,message):

        ''' this method is a placehoder for a possible future exception handler ... notes are below for usage in-line in code '''


        '''        This plain 'except:' catches all exceptions, not only system '''

        #    log = self.log
        #    try:
        #        something_that_may_fail()
        #    except:
        #       error = sys.exc_info()[0]
        #       log.error('ERROR: '+error)


        ''' this is an example of EAFP ... easier to ask for forgiveness than permission ... the preferred way ... also allows some cleanup etc '''

        #def display_username(user_id):
        #    try:
        #        db_connection = get_db_connection()
        #    except DatabaseEatenByGrueError:
        #        print('Sorry! Database was eaten by a grue.')
        #    else:
        #        print(db_connection.get_username(user_id))
        #        db_connection.cleanup()

        ''' the following is an example of check first or look before you leap  LBYL ... could still not work in all cases (race conditions possible) '''

        #def print_object(some_object):
        #    # Check if the object is printable...
        #    try:
        #        printable = str(some_object)
        #    except TypeError:
        #        print("unprintable object")
        #    else:
        #    print(printable)


        ''' see http://www.jeffknupp.com/blog/2013/02/06/write-cleaner-python-use-exceptions/ '''


        ''' proper way to open a file to read to make sure that it closes correctly  if there is an error '''

        #with open("myfile.txt") as f:
        #    for line in f:
        #        print line,

