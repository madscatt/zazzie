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
import os,string,locale,time,platform
import numpy
import Gnuplot,Gnuplot.PlotItems, Gnuplot.funcutils


#       INTERPOLATE
#
#       06/08/2008      --     initial coding				: jc
#
#LC      1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
        INTERPOLATE is the module that calculates an approximate data set
	with a defined number of points and grid spacing.  The output of
	this file is used to compare to synthetic profiles generated in
	subsequent modules.

        This module is called from Data Interpolation from the main
        GUI through the graphical_interpolate.py script.


	REFERENCE:

	Numerical Recipes: The Art of Scientific Computing
	Third Edition (2007), 1256 pp.
	Cambridge University Press
	ISBN-10: 0521880688	
'''

def wait(str=None, prompt='Plot will clear in 10 seconds ...\n'):
        '''
        WAIT is the function to prompt the user to clear a plot on a screen
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
		time.sleep(2)

def print_failure(message,txtOutput):

        txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
        txtOutput.put(">>>> RUN FAILURE <<<<\n")
        txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
        txtOutput.put(message)

        return

def readfile(filename,io,ioe):
    '''
    READFILE is the function to read NCNR SANS data files
    '''

    fake_error = False
    error_magnitude = 0.1
           
    data_file=open(filename,'r').readlines()
#    data_file=open(filename,'r').read().splitlines()

    ### following line replaces control M (\r) with \n

    data_file=[x.replace("\r","\n") for x in data_file]

    x=[] ; y=[] ; z=[] ; nval=1
    x.append(0.0) ; y.append(io); z.append(ioe)	
	
    for line in data_file:
        this_line = string.split(line)
        try:
            qval=locale.atof(this_line[0])
            ival=locale.atof(this_line[1])
            #ival=abs(locale.atof(this_line[1]))
			
            try:
                error_value = locale.atof(this_line[2])
            except:
                error_value = error_magnitude * ival
                fake_error = True
            #if(error_value == 0.0):
            #    error_value = error_value + 1E-5

            x.append(qval)
            y.append(ival)
            z.append(error_value)
            nval=nval+1

        except:
            pass

    if fake_error:
        print 'Error values in I(q) were set to be '+str(100*error_magnitude)+'% of I(0) values since appropriate values were not found'

    print 'x[0] = ',x[0],'; y[0] = ',y[0],'; z[0] = ',z[0]
    print 'x[-1] = ',x[-1],'; y[-1] = ',y[-1],'; z[-1] = ',z[-1]

    print 'len(x) = ',len(x)
    print 'len(y) = ',len(y)
    print 'len(z) = ',len(z)

    return nval,x,y,z

def spline(x,y,nval,yp1,ypn):
        '''
        SPLINE is the function to calculate an approximate data set
        '''

	u=numpy.zeros(nval)
	y2=numpy.zeros(nval)
	maxval=0.99E30 ; nmax=500
	if(yp1>maxval):
		y2[0]=0.0 ; u[0]=0.0
	else:
		y2[0]=-0.5
		u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1)
	for i in range(1,nval-2):
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1])
		p=sig*y2[i-1]+2.0
		y2[i]=(sig-1.0)/p
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
#		u[i]=(6.0*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]))/(x[i+1]-x[i-1])-sig*u[i-1])/p

	if(ypn>maxval):
		qn=0.0 ; un=0.0
	else:
		qn=0.5
		un=(3.0/(x[nval-1]-x[nval-2]))*(ypn-(y[nval-1]-y[nval-2])/(x[nval-1]-x[nval-2]))

	y2[nval-1]=(un-qn*u[nval-2])/(qn*y2[nval-2]+1.0)

	for k in range(nval-2,-1,-1):
		y2[k]=y2[k]*y2[k+1]+u[k]

	return y2

def splint(xa,ya,y2a,nval,x):

	klo=0 ; khi=nval-1 
	while(khi-klo>1):
		if((khi-klo)>1.0):
			k=int((khi+klo)/2)
			if(xa[k]>x): 
				khi=k
			else:
				klo=k
	h=xa[khi]-xa[klo]
	if(h==0.0):
		print 'ERROR: BAD xa INPUT TO ROUTINE SPLIT'
	a=(xa[khi]-x)/h
	b=(x-xa[klo])/h
	ny=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0
	
	return ny

def unpack_variables(variables):

	runname 	= variables['runname'][0]
	expdata		= variables['expdata'][0]
	ofile		= variables['ofile'][0]
	io		= variables['io'][0]    
	ioe		= variables['ioe'][0]    
	dq		= variables['dq'][0]      
	maxpoints	= variables['maxpoints'][0]
	plotflag	= variables['plotflag'][0]

	return runname,expdata,ofile,io,ioe,dq,maxpoints,plotflag

def interpolate(variables,txtOutput):
        '''
        INTERPOLATE is the function to read in variables from GUI input and 
       	calculate an approximate data set to be used in subsequent modeling
	steps. 

        INPUT:  variable descriptions:

		runname:		project name
                expdata:                input NCNR data file (*.sub)
                io:                	I(0) 
                ioe:                	Error in I(0) 
                dq:                  	new delta q 
                maxpoints:             number of new points

        OUTPUT:

                file is stored in "runname"/data_interpolation directory
                
		ofile:                  output filename 

        '''

	runname,expdata,ofile,io,ioe,dq,maxpoints,plotflag=unpack_variables(variables)

        interpath=runname+'/data_interpolation/'
        direxist=os.path.exists(interpath)
        if(direxist==0):
                os.system('mkdir -p '+interpath)

        print 'runname = ',runname

        #ttxt=time.ctime()
        ttxt=time.asctime( time.gmtime( time.time() ) ) 
        st=''.join(['=' for x in xrange(60)])

        txtOutput.put("\n%s \n" %(st))
        txtOutput.put("DATA FROM RUN: %s \n\n" %(ttxt))
        txtOutput.put("Input file used : %s\n" % (expdata))

#	print 'going to attempt to read data'

	nval,x,y,z=readfile(expdata,io,ioe)

#	print 'back from reading data'

	odata=[]
	flag=0; cut=[] 
	for i in range(nval):
		odata.append([x[i],y[i],z[i]])
		if(y[i]/z[i]<2.0 and flag==0):
			cut.append([x[i-1],y[i-1],io])
			cutval=x[i-1]
			flag=1
		elif((i == nval-1) and flag==0):
			cut.append([x[i-1],y[i-1],io])
			cutval=x[i-1]
			flag=1

	yp1=1.0
	ypn=1.0
	y2=numpy.zeros(nval)
	z2=numpy.zeros(nval)

#	print '> calculating splines'

	try:
		y2=spline(x,y,nval,yp1,ypn)
		z2=spline(x,z,nval,yp1,ypn)
	except:
		message='Failed to interpolate data: is data file corrupt?'
		message+='\nstopping here\n\n'
		print_failure(message,txtOutput)
		return

#	print '> back from splines'
	
	outfile2=open(interpath+ofile,'w')
	outfile3=open(interpath+'stn_'+ofile,'w')
       
	if(plotflag == 1): 
		graph = Gnuplot.Gnuplot(debug=1)
        	graph.clear()
        	graph('set title "Interpolation Results"')
        	graph.xlabel('Q (1/A)')
        	graph.ylabel('I(Q)')
		graph('set logscale y')
		
	io_tally=[]
	outfile2.write('%f\t%f\t%f\n' % (0.0,io,ioe))
	outfile3.write('%f\t%f\t%f\n' % (0.0,io,ioe))
	io_tally.append([0.0,io,ioe])
	ux=0.00
	print 'signal to noise cutoff value = ',cutval
	for i in range(maxpoints-1):
		ux=ux+dq
		ny=splint(x,y,y2,nval,ux)
		nz=splint(x,z,z2,nval,ux)
		io_tally.append([ux,ny,nz])
		outfile2.write('%f\t%f\t%f\n' % (ux,ny,nz))
		if(ux<=cutval):
			outfile3.write('%f\t%f\t%f\n' % (ux,ny,nz))
				
	outfile2.close()
	outfile3.close()

	fraction_done=1
        report_string='STATUS\t'+str(fraction_done)
        txtOutput.put(report_string)

        print 'Interpolated data were written to %s\n' % ('./'+interpath+ofile)
        print 'Interpolated data with S/N > 2 were written to %s\n' % ('./'+interpath+'stn_'+ofile)
        print '\ndelta q = %f\t : number of q-points = %i\t : q-range: q = 0 to %f\n' % (dq,maxpoints,(maxpoints-1)*dq)
        txtOutput.put("\nInterpolated data were written to %s\n" % ('./'+interpath+ofile))
        txtOutput.put("\nInterpolated data with S/N > 2 were written to %s\n\n" % ('./'+interpath+'stn_'+ofile))
        txtOutput.put("\ndelta q = %f (1/A)\n\nnumber of q-points = %i\n\nq-range: 0 to %f (1/A)\n" % (dq,maxpoints,(maxpoints-1)*dq))
        txtOutput.put("\n%s \n" %(st))

	if(plotflag == 1): 
		graph.plot(Gnuplot.Data(odata,using='1:2 w p ps 4',title='Original Data'),Gnuplot.Data(io_tally,using='1:2 w lp ps 2',title='Interpolated Data'),Gnuplot.Data(cut,title='[I(Q)/(std.dev. I(Q))] < 2',using='1:2:3 w yerrorbars'))

	time.sleep(2)

	if(plotflag == 1): 
		wait('\n')

	print 'INTERPOLATE IS DONE'
	
	return()



