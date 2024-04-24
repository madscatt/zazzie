import os
import time
import sassie.analyze.hullradsas.HullRadSASV2 as hullradsasv2
import sassie.util.module_utilities as module_utilities
import sassie.util.sasconfig as sasconfig

if sasconfig.__level__ == "DEBUG":
    DEBUG = True

app = 'hullradsas'

class module_variables():

    def __init__(self, parent=None):
        self.app = app

class hullradsas_variables():

    def __init__(self, parent=None):
        pass

class hullradsas():

    def __init__(self, parent=None):
        pass

    def main(self, input_variables, txtOutput):

        self.module_variables = module_variables()

        self.hullradsas_variables = hullradsas_variables()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.run_utils.general_setup(self)

        self.initialization()

        self.run_hullradsas()

        self.epilogue()

        return

    def unpack_variables(self, variables):
        '''
        method to extract variables into system wise class instance
        '''

        log = self.log
        mvars = self.module_variables
        log.debug('in unpack_variables')

        mvars.run_name = variables['run_name'][0]
        mvars.pdbfile = variables['pdbfile'][0]
        mvars.ofile = variables['ofile'][0]

        log.debug(vars(mvars))

        return

    def initialization(self):
        '''
        method to prepare for hullradsas
        '''

#mvars:    run_name, pdbfile, ofile

        log = self.log
        log.debug('in initialization')
        pgui = self.run_utils.print_gui
        mvars = self.module_variables
        hvars = self.hullradsas_variables

        hvars.interpath = mvars.run_name + '/hullradsas/'
        direxist = os.path.exists(hvars.interpath)
        if(direxist == 0):
            os.system('mkdir -p ' + hvars.interpath)

        # ttxt=time.ctime()
        ttxt = time.asctime(time.gmtime(time.time()))
        st = ''.join(['=' for x in range(60)])

        pgui("\n%s \n" % (st))
        pgui("DATA FROM RUN: %s \n\n" % (ttxt))

        return

    def epilogue(self):
        '''
        method to print out results and to move results
        to appropriate places.
        '''

        log = self.log
        mvars = self.module_variables
        hvars = self.hullradsas_variables
        pgui = self.run_utils.print_gui

        log.debug('in epilogue')
        log.debug('HULLRADSAS IS DONE(LOG:EPILOGUE)')


    def run_hullradsas(self):

        log = self.log
        mvars = self.module_variables
        pgui = self.run_utils.print_gui
        log.debug('in run_hullradsas')

        hullradsasv2.main(mvars.pdbfile, mvars.ofile)

        fraction_done = 1
        report_string = 'STATUS\t' + str(fraction_done)
        pgui(report_string)

        self.run_utils.clean_up(log)

        pgui('HULLRADSAS IS DONE')

        pgui("%s \n" % ('=' * 60))
        time.sleep(1.0)
        log.debug('HULLRADSAS IS DONE(LOG)')

        return
