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
import string,os,locale,sys,struct,array,time,platform

try:
    import Gnuplot,Gnuplot.PlotItems, Gnuplot.funcutils
except:
    pass

import numpy
import sassie.util.module_utilities as module_utilities
import sassie.util.sasconfig as sasconfig
import sassie.util.folder_management as folder_management

#       EM_TO_SAS
#
#    12/28/06    --    initial coding              :    jc
#    03/13/09    --    adapted subroutine to SASSIE GUI  :    jc
#    03/31/09    --    added MRC file type input      :    jc
#
#
#LC      1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
        EM_TO_SAS is the module that reads a three-dimensional file
        (either Gaussian Cube or a binary MRC electron density file) and 
        then represent each voxel with an occupation > threshold
        as a scattering center represented as beads (C-atoms), 
        that is then written as PDB/xyz files
        and the SAS profile is calculated using the "scat" program

'''
if sasconfig.__level__ == "DEBUG":
    DEBUG = True

app = 'em_to_sas'

class module_variables():

    def __init__(self, parent=None):
        self.app = app

class em_to_sas_input_variables():

    def __init__(self, parent=None):
        pass


class em_to_sas():

    def __init__(self, parent=None):
        pass

    def main(self, input_variables, txtOutput):

        self.mvars = module_variables()

        self.avars = em_to_sas_input_variables()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.run_utils.general_setup(self)

        self.initialization() 

        self.em()

        self.epilogue()

        return


#    pgui performs this function
#    def print_failure(message,txtOutput):
#
#            txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
#            txtOutput.put(">>>> RUN FAILURE <<<<\n")
#            txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
#            txtOutput.put(message)
#
#            return

    def unpack_variables(self,variables):

        '''
        method to extract variables into system wise class instance
        '''
        log = self.log
        mvars = self.mvars
        log.debug('in unpack_variables')    

        mvars.runname         = variables['runname'][0]
        mvars.emfiletype      = variables['emfiletype'][0]
        mvars.inputpath       = variables['inputpath'][0]
        mvars.emdensityfile   = variables['emdensityfile'][0]
        mvars.pdbfile         = variables['pdbfile'][0]
        mvars.threshold       = variables['threshold'][0]
        mvars.sasfile        = variables['sasfile'][0]
        mvars.npoints         = variables['npoints'][0]
        mvars.qmax            = variables['qmax'][0]
        mvars.plotflag         = variables['plotflag'][0]

#mvars: runname,emfiletype,inputpath,emdensityfile,pdbfile,threshold,sasfile,npoints,qmax,plotflag
#avars: inp, empath, prfile

    def read_mrcfile(self):
        '''
            READ_MRCFILE is the function to read a 3D binary MRC file
            mrc file header parameters are named per the IMOD convention http://bio3d.colorado.edu/imod/doc/mrc_format.txt
        '''

        log = self.log
        mvars = self.mvars
        avars = self.avars
        pgui = self.run_utils.print_gui
        log.debug('in read_mrcfile')

#       no longer checked; mrc file is read the same way using both 32-bit and 64-bit Python
#       bit = platform.architecture()[0]

        ang2au=1.0/0.5291772108
        au2ang=1.0/ang2au
        outfile=open(mvars.inputpath+'/'+mvars.pdbfile,'w')
        outfile2=open(mvars.inputpath+'/'+mvars.pdbfile+'.xyz','w')

        input=open(mvars.emdensityfile,'rb').read()

#       read the header

        start,stop=0,struct.calcsize('iiiiiiiiiiffffffiiifffii')    
        nc,nr,ns,mode,ncstart,nrstart,nsstart,mx,my,mz,xlen,ylen,zlen,alpha,beta,gamma,mapc,mapr,maps,amin,amax,amean,ispg,nsymbt=struct.unpack('iiiiiiiiiiffffffiiifffii',input[start:stop])    

        pgui('nc, nr, ns: %5i%5i%5i' % (nc,nr,ns))
        pgui('mode/type: %i'% (mode))
        pgui('ncstart, nrstart, nsstart: %5i%5i%5i'% (ncstart,nrstart,nsstart))
        pgui('xlen, ylen, zlen: %8.3f%8.3f%8.3f' % (xlen, ylen, zlen))
        pgui('unit cell size (mx, my, mz):  %5i%5i%5i' % (mx,my,mz))
        pgui('cell angles: %5.1f%5.1f%5.1f' % (alpha, beta, gamma))
        pgui('mapc, mapr, maps: %5i%5i%5i' % (mapc, mapr, maps))        #used to compute axes permutation (Chimera)
        pgui('maximum intensity: %f' % (amax))
        pgui('minimum intensity: %f' % (amin))
        pgui('amean: %f' % (amean))
        pgui('ispg: %i' % (ispg))
        pgui('nsymbt: %i' % (nsymbt))

#determine data step and angstroms per pixel resolution from header parameters

        if mx > 0 and my > 0 and mz > 0 and xlen > 0 and ylen > 0 and zlen > 0:
            xdatastep = mx/xlen
            ydatastep = my/ylen
            zdatastep = mz/zlen
            xangpix = xlen/mx
            yangpix = ylen/my
            zangpix = zlen/mz
            pgui('data_step (from header parameters): %8.3f%8.3f%8.3f' % (xdatastep, ydatastep, zdatastep))
            pgui('angstroms/pixel resolution (from header parameters): %8.3f%8.3f%8.3f' % (xangpix,yangpix,zangpix))
        else:
            xdatastep = 1.0
            ydatastep = 1.0
            zdatastep = 1.0        
            xangpix = 1.0
            yangpix = 1.0
            zangpix = 1.0
            pgui('data_step (default): %8.3f%8.3f%8.3f' % (xdatastep, ydatastep, zdatastep))
            pgui('angstroms/pixel resolution (default): %8.3f%8.3f%8.3f' % (xangpix,yangpix,zangpix))

#determine correct x,y,z origin for possible future use

        type_string = struct.unpack('ssss',input[208:212])    
        pgui('type_string = %s' % (str(type_string)))
        if type_string == ('M', 'A', 'P', ' '):
            filetype = 'mrc2000'        #newest mrc format
            xyz_origin = struct.unpack('fff',input[196:208])
            pgui('xyz_origin: %s' % (str(xyz_origin)))
        else:
            filetype = 'mrc'
            xyz_origin = (ncstart,nrstart,nsstart)

        xorigin = xyz_origin[0]
        yorigin = xyz_origin[1]
        zorigin = xyz_origin[2]
        
        pgui('file type: %s' % (filetype))
        pgui('data_origin: %5i%5i%5i' % (xorigin, yorigin, zorigin))
    
#   read the data

        start=1024
        fsize=struct.calcsize('f')
        stop=start+fsize
        pgui('fsize: %i' % (fsize))
        pgui('start: %i' % (start))
        pgui('stop: %i' % (stop))

        nxgp=mx ; nygp=my ; nzgp=mz
        zmin=-mz/2.0 ; ymin=-my/2.0 ; xmin=-mx/2.0
        xgs=(2.0*abs(xmin)/nxgp)*xangpix ; xgs2=xgs/2.0
        ygs=(2.0*abs(ymin)/nygp)*yangpix; ygs2=ygs/2.0
        zgs=(2.0*abs(zmin)/nzgp)*zangpix ; zgs2=zgs/2.0
        zmina=zmin*zangpix ; ymina=ymin*yangpix ; xmina=xmin*xangpix    
        ii=0
        rbin=numpy.zeros((nxgp*nygp*nzgp,3),numpy.float32)
        for i in range(mz):
            thisz=zmina+zgs2+i*zgs
            for j in range(my):
                thisy=ymina+ygs2+j*ygs
                for k in range(mx):
                    thisx=xmina+xgs2+k*xgs
                    val=struct.unpack('f',input[start:stop])
                    start=stop
                    stop=stop+fsize
                    if(val[0]>mvars.threshold):
                        rbin[ii][0]=thisx    
                        rbin[ii][1]=thisy    
                        rbin[ii][2]=thisz    
                        ii=ii+1
    
        numi=ii
    
        st1="ATOM  "
        st2="  C   DUM"
        st3="    "
        st4="  0.00  0.00      DUM"
#ATOM      2  HT1 LYS     1     -13.587 -15.185  12.469  0.32  0.60      7RSA
        outfile.write("REMARK\n")
    
#For scat, there is no need to have more than about 10000 points in the volume.  So, the number of points is reduced to 10000 for testing purposes while using scat to calculate the scattering.  

        for i in range(numi):

            if(i>9999):
                resid = 9999
            else:
                resid = i+1

            if(i>99999):
                index = 99999
            else:
                index = i+1

            outfile.write('%6s%5i%s %5i%s%8.3f%8.3f%8.3f%s\n' % (st1,index,st2,resid,st3,rbin[i][0],rbin[i][1],rbin[i][2],st4))
            outfile2.write('%f\t%f\t%f\t%f\n' % (rbin[i][0],rbin[i][1],rbin[i][2],1.0))
        outfile.close()    
        outfile2.close()    
    
        return

    def wait(self,str=None, prompt='Plot will clear in 10 seconds ...\n'):
        '''
            WAIT is the function to wait for user input to proceed.
        '''
    
        if str is not None:
            print str

        try:
            if(platform.system() == "Linux"):
                import curses
                stdscr = curses.initscr()
                stdscr.addstr('press a key to continue')
                c = stdscr.getch()
                curses.endwin()
        except:
            time.sleep(1)

    def read_cubefile(self):
        '''
            READ_CUBEFILE is the function to read a 3D binary Gaussian cube file
        '''
    
        log = self.log
        mvars = self.mvars
        avars = self.avars
        pgui = self.run_utils.print_gui
        log.debug('in read_cubefile')

        ang2au=1.0/0.5291772108
        au2ang=1.0/ang2au
    
        outfile=open(mvars.inputpath+'/'+mvars.pdbfile,'w')
        outfile2=open(mvars.inputpath+'/'+mvars.pdbfile+'.xyz','w')
        infile=open(mvars.emdensityfile,'r').readlines()

        line0=infile[0]        #calpha_JC CUBE FILE
        line1=infile[1]        #OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z
        line2=infile[2]        #431     -236.215766     -236.215766     -236.215766
        line3=infile[3]        #41      11.338357       0.000000        0.000000
        line4=infile[4]        #41      0.000000        11.338357       0.000000
        line5=infile[5]        #41      0.000000        0.000000        11.338357    
        lin2=string.split(line2)
        lin3=string.split(line3)
        lin4=string.split(line4)
        lin5=string.split(line5)

        natom=locale.atoi(lin2[0])
        xmin=au2ang*locale.atof(lin2[1])    
        ymin=au2ang*locale.atof(lin2[2])    
        zmin=au2ang*locale.atof(lin2[3])    

        nxgp=locale.atoi(lin3[0])
        nygp=locale.atoi(lin4[0])
        nzgp=locale.atoi(lin5[0])

        pgui('natom: %i' % (natom))
        pgui('xmin, ymin,zmin: %8.3f%8.3f%8.3f' % (xmin,ymin,zmin))
        pgui('nxgp, nygp, nzgp:  %8.3f%8.3f%8.3f' % (nxgp,nygp,nzgp)) 

        cbin=numpy.zeros((nxgp,nygp,nzgp),numpy.float32)
        rbin=numpy.zeros((nxgp*nygp*nzgp,3),numpy.float32)
    
        xgs=(2.0*abs(xmin)/nxgp) ; xgs2=xgs/2.0
        ygs=(2.0*abs(ymin)/nygp) ; ygs2=ygs/2.0
        zgs=(2.0*abs(zmin)/nzgp) ; zgs2=zgs/2.0

        pgui('xgs, ygs, zgs:  %8.3f%8.3f%8.3f' % (xgs,ygs,zgs))

        skip=natom+6 
        line=infile[skip] ; thisline=skip
        lin=string.split(line) 
        ii=0    
        for i in range(nxgp):
            thisx=xmin+xgs2+i*xgs
            for j in range(nygp):
                thisy=ymin+ygs2+j*ygs
                loc=0
                for k in range(nzgp):
                    thisz=zmin+ygs2+k*zgs
                    l=len(lin)
                    cbin=locale.atof(lin[loc])
                    stop=0    
                    if(i==nxgp-1 and j==nygp-1 and k==nzgp-1):
                        stop=1
                    if(((k%6)==5 or k==nzgp-1) and stop==0) :
                        thisline=thisline+1
                        line=infile[thisline]    
                        lin=string.split(line)
                        loc=0
                    else:
                        loc=loc+1

                    if(cbin>mvars.threshold):
                        rbin[ii][0]=thisx    
                        rbin[ii][1]=thisy    
                        rbin[ii][2]=thisz    
                        ii=ii+1

        numi=ii

        st1="ATOM  "
        st2="  C   DUM"
        st3="    "
        st4="  0.00  0.00      DUM"
#ATOM      2  HT1 LYS     1     -13.587 -15.185  12.469  0.32  0.60      7RSA
        outfile.write("REMARK\n")
        for i in range(numi):
            outfile.write('%6s%5i%s %5i%s%8.3f%8.3f%8.3f%s\n' % (st1,i+1,st2,i+1,st3,rbin[i][0],rbin[i][1],rbin[i][2],st4))
            outfile2.write('%f\t%f\t%f\t%f\n' % (rbin[i][0],rbin[i][1],rbin[i][2],1.0))
        outfile.close()    
        outfile2.close()    
    
        return

    def write_inputfile(self):
        '''
            WRITE_INPUTFILE is the function to write the input file for scat
        '''
    
        log = self.log
        mvars = self.mvars
        avars = self.avars
        log.debug('in write_inputfile')


        outfile=open(avars.inp,'w')
        outfile.write('%i\n%s\n%s\n%i\n%f\n' % (1,mvars.pdbfile+'.xyz',mvars.sasfile,mvars.npoints,mvars.qmax))
        outfile.close()    

        return

    def run_scat(self):
        '''
            RUN_SCAT is the function to run the scat program
        '''

        log = self.log
        mvars = self.mvars
        avars = self.avars
        log.debug('in run_scat')


        bin_path = sasconfig.__bin_path__
    
        currcommand=bin_path+'/scat.exe < '+avars.inp
        log.debug('currcommand: %s' %(currcommand))
        os.system(currcommand)

        return
 
    def read_dat(self,ifile):
        '''
            READ_DAT is the function to read data from the input I(q) or P(r) file
        '''

        infile=open(ifile,'r').readlines()
        nl=len(infile)
        data=[]
        for i in range(nl):
            lin=string.split(infile[i])
            x=locale.atof(lin[0]) ; y=locale.atof(lin[1]) ; z=locale.atof(lin[2])
            data.append([x,y,z])

        return data

    def initialization(self):
        '''
        method to prepare for em
        '''

        log = self.log
        log.debug('in initialization')
        pgui = self.run_utils.print_gui

        mvars = self.mvars
        avars = self.avars

        avars.empath=mvars.runname+'/em_to_sas'
        direxist=os.path.exists(avars.empath)
        if(direxist==0):
            try:
                result=os.system('mkdir -p '+avars.empath)
            except:
                message='can not create project directory: '+avars.empath
                message+='\nstopping here\n'
                pgui(message)
            if(result!=0):
                message='can not create project directory: '+avars.empath
                message+='\nstopping here\n'
                pgui(message)

        return
        

    def em(self):

        '''
            EM is the function to read in three-dimensional voxel data and
            calculate the SAS profile through the binary program scat.

            INPUT:  variable descriptions:

                runname:        project name
                emfiletype:     em file type (0=cube,1=mrc)
                inputpath:      input file path
                emdensityfile:  em filename
                threshold:      treshold cutoff
                npoints:        number of points in sas calculation
                qmax:           q-max in sas calculation

            OUTPUT:

                pdbfile:        output filename (pdb)
                
                files stored in ./"runname"/em_to_sas/ directory:

                sasfile*.sub:  output sas profile
                sasfile*.pr:   output p(r)
                dum.inp:       input file written for scat
                sasfile*.pdb:  pdb file of coordinates used for scat
                sasfile*.xyz:  xyz file of coordinates used for scat
        '''


        log = self.log
        log.debug('in em')
        pgui = self.run_utils.print_gui

        mvars = self.mvars
        avars = self.avars

        log.debug(vars(mvars))
        log.debug(vars(avars))


        pgui('runname = %s' % (mvars.runname))
        
        #ttxt=time.ctime()
        ttxt=time.asctime( time.gmtime( time.time() ) ) 
        st=''.join(['=' for x in xrange(60)])

        pgui("\n%s \n" %(st))
        pgui("DATA FROM RUN: %s \n\n" %(ttxt))

        avars.inp='dum.inp'
        if(mvars.emfiletype==0):
            self.read_cubefile()
        elif(mvars.emfiletype==1):
            self.read_mrcfile()
#this error message should not occur since emfiletype is check in the em_to_sas filter
        else:
            pgui('wrong file type entered')
            pgui('0==Gaussian Cube file')
            pgui('1==MRC file')
            pgui('stopping now')
            message='wrong file type entered: '+mvars.emfiletype
            message+='\nstopping here\n'
            pgui(message)
            sys.exit()

        fraction_done = 0.5 
        report_string='STATUS\t'+str(fraction_done)
        pgui(report_string)

        self.write_inputfile()
        self.run_scat()
    
        fraction_done = 1.0 
        report_string='STATUS\t'+str(fraction_done)
        pgui(report_string)


        if(mvars.plotflag == 1):
            graph = Gnuplot.Gnuplot(debug=1)
            graph.clear()
            graph('set title "SAS Profile"')
            graph.xlabel('Q (1/A)')
            graph.ylabel('I(Q)')
            graph('set logscale x')
            graph('set logscale y')

            iqdat=self.read_dat(mvars.sasfile)

            graph.plot(Gnuplot.Data(iqdat,using='1:2 w lp ps 4',title='Calculated SAS Profile'))

        for i in range(len(mvars.sasfile)):
            char=mvars.sasfile[i]
            if(char=='.'):
                pos=i
        avars.prfile=mvars.sasfile[0:pos]+'.pr'
        log.debug('prfile: %s' % (avars.prfile))

        prdat=self.read_dat(avars.prfile)

        pgui("Data stored in directory: %s\n\n" % ('./'+avars.empath))
        pgui("\n%s \n" %(st))
        time.sleep(1)

        if(mvars.plotflag == 1):
            graph2 = Gnuplot.Gnuplot(debug=1)
            graph2.clear()
            graph2('set title "P(r) Profile"')
            graph2.xlabel('r (Angstrom)')
            graph2.ylabel('P(r)')

            graph2.plot(Gnuplot.Data(prdat,using='1:2 w lp ps 4',title='Calculated P(r) Profile'))
            
        if(mvars.plotflag == 1):
            self.wait('\n')

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

        log.debug('pdbfile: %s' %(mvars.pdbfile))
        mvst='mv -f '+avars.prfile+' '+mvars.pdbfile+' '+avars.inp+' '+mvars.sasfile+' '+mvars.pdbfile+'.xyz '+avars.empath
        log.debug('mvst: %s' %(mvst))
        os.system(mvst)

        self.run_utils.clean_up(log)

        pgui("%s \n" % ('=' * 60))
        pgui('EM TO SAS IS DONE')
        time.sleep(1.5)

        return        









