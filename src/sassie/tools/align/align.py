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
import os
import sys
import string
import locale
import time
import numpy
import sasmol.sasmol as sasmol
import sassie.util.module_utilities as module_utilities
import sassie.util.sasconfig as sasconfig

'''
    ALIGN is the module that overlaps molecules from a dcd/pdb file
    onto another molecule over a given basis.  The two molecule types
    do not need to be the same but the basis atoms and the number of
    basis atoms used for the overlap do need to be identical.

    REFERENCE:

    W. Kabsch
    Acta Crystallog. sect. A  32  922-923  (1976)

    W. Kabsch
    Acta Crystallog. sect. A  34  827-828  (1978)    
'''


if sasconfig.__level__ == "DEBUG":
    DEBUG = True

app = 'align'

class module_variables():

    def __init__(self, parent=None):
        self.app = app

class align_input_variables():

    def __init__(self, parent=None):
        pass

class align():

    def __init__(self, parent=None):
        pass

    def main(self, input_variables, txtOutput):

        self.mvars = module_variables()

        self.avars = align_input_variables()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.run_utils.general_setup(self)

        self.initialization()

        self.align(input_variables, txtOutput)

        self.epilogue()

        return

    def unpack_variables(self, variables):
        '''
        method to extract variables into system wise class instance
        '''

        log = self.log
        mvars = self.mvars
        log.debug('in unpack_variables')

        mvars.runname = variables['runname'][0]
        mvars.path = variables['path'][0]
        mvars.infile = variables['infile'][0]
        mvars.pdbmol1 = variables['pdbmol1'][0]
        mvars.pdbmol2 = variables['pdbmol2'][0]
        mvars.ofile = variables['ofile'][0]
        mvars.basis1 = variables['basis1'][0]
        mvars.basis2 = variables['basis2'][0]
        mvars.lowres1 = variables['lowres1'][0]
        mvars.lowres2 = variables['lowres2'][0]
        mvars.highres1 = variables['highres1'][0]
        mvars.highres2 = variables['highres2'][0]
        mvars.ebasis1 = variables['ebasis1'][0]
        mvars.ebasis2 = variables['ebasis2'][0]

        log.debug(vars(mvars))

        return

#    pgui performs this function 
#    def print_failure(self,message,txtOutput):
#
#        txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
#        txtOutput.put(">>>> RUN FAILURE <<<<\n")
#        txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
#        txtOutput.put(message)
#
#        return


    def write_frame_to_file(self,frame):        

        mvars = self.mvars
        avars = self.avars

        if(avars.outtype == 'dcd'):
            if(avars.intype == 'dcd'):
                avars.m2.write_dcd_step(avars.dcdoutfile, 0, frame)
            elif(avars.intype == 'pdb'):
                if(avars.nf==1):
                    avars.m2.write_dcd_step(avars.dcdoutfile,0,frame)
                else:
                    avars.m2.write_dcd_step(avars.dcdoutfile,frame,frame)

        elif(avars.outtype == 'pdb'):
            pdbfilestring = os.path.join(avars.alignpath,mvars.ofile)
            if(avars.nf == 1):
                avars.m2.write_pdb(pdbfilestring,0,'w')
            elif(avars.nf > 1 and avars.intype == 'dcd'):
                avars.m2.write_pdb(pdbfilestring,0,'a')
            elif(avars.nf > 1 and avars.intype == 'pdb'):
                avars.m2.write_pdb(pdbfilestring,frame,'a')


    def initialization(self):
        '''
        method to prepare for alignment
        '''

        log = self.log
        log.debug('in initialization')
        pgui = self.run_utils.print_gui
        mvars = self.mvars
        avars = self.avars


        ''' directory and file preparation '''

        avars.alignpath=os.path.join(mvars.runname,app)
        log.debug('align path: '+avars.alignpath)
        direxist=os.path.exists(avars.alignpath)
        if(direxist == 0):
            try:
                result = os.system('mkdir -p ' + avars.alignpath)
            except:
                message = 'can not create project directory: ' + avars.alignpath
                message += '\nstopping here\n'
                pgui(message)
                sys.exit(1)
            if(result != 0):
                message = 'can not create project directory: ' + avars.alignpath
                message += '\nstopping here\n'
                pgui(message)
                sys.exit(1)

        if(mvars.ofile[-3:] == 'dcd'):
            log.debug('output file is a DCD file')
            avars.outtype = 'dcd'
        elif(mvars.ofile[-3:] == 'pdb'):
            log.debug('output file is a PDB file')
            avars.outtype = 'pdb'

#        this condition is checked in align_filter.py and will stop if there is no file extension, i.e., can't get to here
#        else:
#            avars.outtype = 'dcd'
#            message='output filename '+mvars.ofile+' needs to end in either ".pdb" (1 frame) or ".dcd" (1 or more frames)\n'
#            message+=' :  writing output file as a '+mvars.ofile+'.dcd\n'
#            pgui('\n\n',message,'\n\n')
#            mvars.ofile = mvars.ofile +'.dcd'

        avars.minmaxfile=mvars.ofile+'.minmax'

        avars.m1=sasmol.SasMol(0)
        avars.m2=sasmol.SasMol(1)

        avars.m1.read_pdb(mvars.path+mvars.pdbmol1,check_zero_coor=True)
        avars.m2.read_pdb(mvars.path+mvars.pdbmol2,check_zero_coor=True)

        if(mvars.infile[-3:] == 'dcd'):
            dcdfile = avars.m2.open_dcd_read(mvars.path+mvars.infile)
            avars.nf = dcdfile[2]
            avars.intype = 'dcd'
            log.debug('>> input file is a DCD file')
            log.debug('number of frames: '+avars.nf)
        elif(mvars.infile[-3:] =='pdb'):
            avars.m2.read_pdb(mvars.path+mvars.infile)
            avars.nf = avars.m2.number_of_frames()
            avars.intype = 'pdb'
            log.debug('input file is a PDB file')
            log.debug('number of frames: '+str(avars.nf))			


#        this condition is checked in align_filter.py
#        else:
#            message='input filename is a PDB or DCD file but it must end with ".pdb" or ".dcd" '
#            message+=' :  stopping here'
#            pgui(message)
#            sys.exit(1)


        avars.dcdoutfile = 'None'	
        if(avars.outtype == 'dcd'):
            dcdfilestring = os.path.join(avars.alignpath,mvars.ofile)
            log.debug('dcdoutfile: '+dcdfilestring)
            avars.dcdoutfile = avars.m2.open_dcd_write(dcdfilestring)
            pgui("Total number of frames = %d\n\n" % (avars.nf))


        mass1=avars.m1.mass()
        mass2=avars.m2.mass()

        name1=avars.m1.name()
        name2=avars.m2.name()

        basis_filter_1_a = '(name[i] == "'+mvars.basis1+'" and (resid[i] >= '+str(mvars.lowres1)+' and resid[i] <= '+str(mvars.highres1)+'))'
        basis_filter_2_a = '(name[i] == "'+mvars.basis2+'" and (resid[i] >= '+str(mvars.lowres2)+' and resid[i] <= '+str(mvars.highres2)+'))'

        if(mvars.ebasis1 == "None"):
            basis_filter_1 = basis_filter_1_a
        else:
            basis_filter_1 = ebasis1 + ' and '+basis_filter_1_a

        if(mvars.ebasis2 == "None"):
            basis_filter_2 = basis_filter_2_a	
        else:
            basis_filter_2 = ebasis2 + ' and '+basis_filter_2_a

        error,avars.mask1 = avars.m1.get_subset_mask(basis_filter_1)	
        error,avars.mask2 = avars.m2.get_subset_mask(basis_filter_2)	

        avars.sub_m1 = sasmol.SasMol(2)
        error = avars.m1.copy_molecule_using_mask(avars.sub_m1,avars.mask1,0)

        avars.sub_m2 = sasmol.SasMol(3)
        error = avars.m2.copy_molecule_using_mask(avars.sub_m2,avars.mask2,0)

        avars.com_sub_m1=avars.sub_m1.calccom(0)
        avars.sub_m1.center(0)
        avars.coor_sub_m1=avars.sub_m1.coor()[0]

        return


    def align(self,variables,txtOutput):
        '''
        ALIGN is the function to read in variables from GUI input and 
        overlap the molecules in a dcd/pdb file onto the coordinates of
        a reference pdb structure over a given basis.

        INPUT:  variable descriptions: 

        runname: 	            project name
        path:                   input/output filepath
        pdbmol1:                reference pdb (mol 1)
        pdbmol2:                input pdb file (mol 2)
        infile:                 input (pdb or dcd) filename (mol 2)
        basis1:                 basis for molecule 1
        basis2:                 basis for molecule 2
        lowres1:                low residue for overlap molecule 1
        highres1:               high residue for overlap molecule 1
        lowres2:                low residue for overlap molecule 2
        highres2:               high residue for overlap molecule 2
        ebasis1:		            extra basis statement molecule 1 
        ebasis2:		            extra basis statement molecule 2 

        OUTPUT: files stored in "runname"/align directory:

        ofile:			        output filename
        ofile*.minmax:		    text file with min & max dimensions
        '''
        log = self.log
        pgui = self.run_utils.print_gui
        write_frame_to_file = self.write_frame_to_file
        log.debug('in align')

        mvars = self.mvars
        avars = self.avars

        minmaxfilestring = os.path.join(avars.alignpath,avars.minmaxfile)
        log.debug('minmaxfile: '+minmaxfilestring)
        mmfile=open(minmaxfilestring,'w')

        ttxt=time.asctime( time.gmtime( time.time() ) ) 
        st=''.join(['=' for x in xrange(60)])
        pgui("\n%s \n" %(st))
        pgui("DATA FROM RUN: %s \n\n" %(ttxt))



        minx=[] ; miny=[] ; minz=[]
        maxx=[] ; maxy=[] ; maxz=[]

        for i in xrange(avars.nf):
            if(avars.intype == 'dcd'):
                avars.m2.read_dcd_step(avars.dcdfile,i)
                avars.m2.center(0)
                minmax = avars.m2.calcminmax_frame(0)
                error,avars.sub_m2.coor = avars.m2.get_coor_using_mask(0,avars.mask2)
                avars.sub_m2.setCoor(avars.sub_m2.coor)
                avars.com_sub_m2 = avars.sub_m2.calccom(0)
                avars.sub_m2.center(0)
                avars.coor_sub_m2 = avars.sub_m2.coor[0]
                avars.m2.align(0,avars.coor_sub_m2,avars.com_sub_m2,avars.coor_sub_m1,avars.com_sub_m1)
            elif(avars.intype == 'pdb'):
                avars.m2.center(i)
                minmax = avars.m2.calcminmax_frame(i)
                error,avars.sub_m2.coor = avars.m2.get_coor_using_mask(i,avars.mask2)
                avars.sub_m2.setCoor(avars.sub_m2.coor)
                avars.com_sub_m2 = avars.sub_m2.calccom(0)
                avars.sub_m2.center(0)
                avars.coor_sub_m2 = avars.sub_m2.coor[0]
                avars.m2.align(i,avars.coor_sub_m2,avars.com_sub_m2,avars.coor_sub_m1,avars.com_sub_m1)
            log.debug(vars(avars))
            

            write_frame_to_file(i+1)

            minx.append(minmax[0][0]) ; miny.append(minmax[0][1]) ; minz.append(minmax[0][2])
            maxx.append(minmax[1][0]) ; maxy.append(minmax[1][1]) ; maxz.append(minmax[1][2])

            ''' display progress '''
            
            if(((i+1)%(float(avars.nf)/10.0)==0 or (avars.nf<10))):
                fraction_done = (float(i+1)/float(avars.nf))
                progress_string='\nCOMPLETED '+str(i+1)+' of '+str(avars.nf)+' : '+str(fraction_done*100.0)+' % done'
                pgui('%s\n\n' % progress_string)
                report_string='STATUS\t'+str(fraction_done)
                log.debug(report_string)

        if(avars.intype == 'dcd'):
            avars.m2.close_dcd_read(avars.dcdfile[0])
        if(avars.outtype == 'dcd'):
            avars.m2.close_dcd_write(avars.dcdoutfile)

        avars.min_x = numpy.min(minx) ; avars.min_y = numpy.min(miny) ; avars.min_z = numpy.min(minz)
        avars.max_x = numpy.max(maxx) ; avars.max_y = numpy.max(maxy) ; avars.max_z = numpy.max(maxz)

        mmfile.write("%s\n" % ("#min_x,min_y,min_z,max_x,max_y,max_z"))
        mmfile.write("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n" % (avars.min_x,avars.min_y,avars.min_z,avars.max_x,avars.max_y,avars.max_z))
        mmfile.close()

        return

    def epilogue(self):
        '''
        method to print out results and to move results
        to appropriate places.
        '''

        log = self.log
        mvars = self.mvars
        avars = self.avars
        pgui = self.run_utils.print_gui

        log.debug('in epilogue')

        pgui("minimum x = %lf\t maximum x = %lf -> range: %lf Angstroms\n" % (avars.min_x,avars.max_x,(avars.max_x-avars.min_x)))
        pgui("minimum y = %lf\t maximum y = %lf -> range: %lf Angstroms\n" % (avars.min_y,avars.max_y,(avars.max_y-avars.min_y)))
        pgui("minimum z = %lf\t maximum z = %lf -> range: %lf Angstroms\n\n" % (avars.min_z,avars.max_z,(avars.max_z-avars.min_z)))

        filestring = os.path.join(avars.alignpath,mvars.ofile)
        minmaxfilestring = os.path.join(avars.alignpath,avars.minmaxfile)
        pgui("\nAligned data (nf=%i) were written to %s\n" % (avars.nf,filestring)) 
        pgui("\nDimension data were written to %s\n" % (minmaxfilestring))

        self.run_utils.clean_up(log)

        pgui("%s \n" % ('=' * 60))
        time.sleep(1.0)

        return

