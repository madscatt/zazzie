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
import os,sys,time,string,locale,random,numpy,time,platform
import Gnuplot,Gnuplot.PlotItems, Gnuplot.funcutils

#       CHI_SQUARE
#
#       12/5/2004       --      initial coding                  :       jc
#
#LC      1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
        CHI_SQUARE is the main module that contains the functions
      	that are used to compare interpolated experimental data
	to the synthetic data sets generated from structures files
	using the modules in CALCULATE.

	Some filtering options based on Rg are provided.

	This module is called from Chi-Square Filter from the main
	GUI through the graphical_chi_square.py script.
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
		time.sleep(1)

def print_failure(message,txtOutput):

        txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
        txtOutput.put(">>>> RUN FAILURE <<<<\n")
        txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
        txtOutput.put(message)

        return

def unpack_variables(variables):

        runname                 = variables['runname'][0]
        saspath                 = variables['saspath'][0]
        sasintfile              = variables['sasintfile'][0]
        io                      = variables['io'][0]
        x2highcut               = variables['x2highcut'][0]
        x2highweightfile        = variables['x2highweightfile'][0]
        x2lowcut                = variables['x2lowcut'][0]
        x2lowweightfile         = variables['x2lowweightfile'][0]
        rghighcut               = variables['rghighcut'][0]
        rghighweightfile        = variables['rghighweightfile'][0]
        rglowcut                = variables['rglowcut'][0]
        rglowweightfile         = variables['rglowweightfile'][0]
        sastype                 = variables['sastype'][0]
        reduced_x2 		= variables['reduced_x2'][0]
        plotflag			= variables['plotflag'][0]

	return runname,saspath,sasintfile,io,x2highcut,x2highweightfile,x2lowcut,x2lowweightfile,rghighcut,rghighweightfile,rglowcut,rglowweightfile,sastype,reduced_x2,plotflag

def chi_square(variables,txtOutput):

	'''
        CHI_SQUARE is the function to read in variables from GUI input and compare
	interpolated experimental data to synthetic data sets generated from
	structure files using the modules in CALCULATE.

	INPUT:	variable descriptions:

                runname:               used for automatic path generation
                saspath:                path to synthetic SAS files
                sasintfile:             name of interpolated data set
                io:                    I(0) value
                x2highcut:             X2 high cutoff value
                x2lowcut:              X2 low cutoff value
                rghighcut:             Rg high cutoff value
                rglowcut:              Rg low cutoff value
                sas_type:               designate scattering calculator (1=xtal2sas, 2=cryson, 3=crysol)
	          reduced_x2:		    set option for using reduced X2 (1==YES, 0==NO)
		
        OUTPUT:
                x2highweightfile:       name of file containing weights for X2 high cutoff
                x2lowweightfile:        name of file containing weights for X2 low cutoff
                rghighweightfile:       name of file containing weights for Rg high cutoff
                rglowweightfile:        name of file containing weights for Rg low cutoff
                txtOutput:              TK handler for output to GUI textbox

                files stored in ~/runname/filter directory:

                bestworstfile.txt:      best and worst sas profiles from single structure
                averagefile.txt:        average profile for the entire ensemble
                x2file.txt:             list of X2 and Rg with each filename
                x2_vs_rg_plot.txt:      Rg and X2 without filename (for plotting)
                sas_spectra_plot.txt:   text file of sas profile plot from Gnuplot
	'''

	runname,saspath,sasintfile,io,x2highcut,x2highweightfile,x2lowcut,x2lowweightfile,rghighcut,rghighweightfile,rglowcut,rglowweightfile,sastype,reduced_x2,plotflag=unpack_variables(variables)

	epsilon=0.00001

	bestworstfile='bestworstfile.txt'
	averagefile='averagefile.txt'
	x2file='x2file.txt'

        filter=runname+'/filter'
        direxist=os.path.exists(filter)
        if(direxist==0):
                os.system('mkdir -p '+filter)
 
        spectra=runname+'/filter/spectra'
        direxist=os.path.exists(spectra)
        if(direxist==0):
                os.system('mkdir -p '+spectra) 
  
	sasfilestring=runname+'_'

	ttxt=time.ctime()

	ginfile=open(sasintfile,'r').readlines()
	ngpoints=len(ginfile)
	lowlin=string.split(ginfile[0])
	low=locale.atof(lowlin[0])
	lowplusonelin=string.split(ginfile[1])
	lowplusone=locale.atof(lowplusonelin[0])
	highlin=string.split(ginfile[-1])
	high=locale.atof(highlin[0])
	
	step=lowplusone-low

	x2cl=[] ; x2ch=[] ; rgcl=[] ; rgch=[]
	runst=runname+'/filter/'
	x2lowfile=open(runst+x2lowweightfile,'w')
	x2highfile=open(runst+x2highweightfile,'w')
	rglowfile=open(runst+rglowweightfile,'w')
	rghighfile=open(runst+rghighweightfile,'w')
	
	gridq=[]
	igoal=[]
	isigma=[]
	for i in xrange(ngpoints):
		lin=string.split(ginfile[i])
		gval=locale.atof(lin[0])
		val=locale.atof(lin[1])
		val2=locale.atof(lin[2])
		gridq.append(gval)
		igoal.append(val)
		isigma.append(val2)

	goal=numpy.array(igoal)
	sigma=numpy.array(isigma)
	sigma2=sigma*sigma
	
	outfile=open(runst+averagefile,'w')
	x2outfile=open(runst+x2file,'w')
	bestworstoutfile=open(runst+bestworstfile,'w')
	
	nfilestring='find '+saspath+' -name "*.log" | grep -c log'
	nfout=os.popen(nfilestring,'r').readlines()
	nf=locale.atoi(nfout[0])
	print "\n# READING Rg VALUES FROM LOG FILES"	
	log=[] 

	for i in xrange(nf):

		mst = str(i+1).zfill(5)  #99999 maximum number of frames   
        	log.append(saspath+'/'+sasfilestring+mst)
			
	ns=nf		# number of trajectory steps
	rgarray=[] ;keep=[] ; allx2=[] 
	for i in xrange(ns):
		inf=open(log[i]+'.log','r').readlines()
		ln=len(inf)
		for k in xrange(ln):
			lin=string.split(inf[k])
			if(len(lin)>0):
				if(sastype==1):
					if(lin[1]=='Rg'):
						trg=lin[3]
						ftrg=locale.atof(trg)
						keep.append(i) 
						rgarray.append(ftrg)
						if(ftrg<rglowcut):
							rgcl.append(1.0)
						else:
							rgcl.append(0.0)
						if(ftrg>rghighcut):
							rgch.append(1.0)
						else:
							rgch.append(0.0)
				if(sastype==2):
					if(lin[0]=='Atomic'):
						trg=lin[3]
						ftrg=locale.atof(trg)
						keep.append(i) 
						rgarray.append(ftrg)
						if(ftrg<rglowcut):
							rgcl.append(1.0)
						else:
							rgcl.append(0.0)
						if(ftrg>rghighcut):
							rgch.append(1.0)
						else:
							rgch.append(0.0)
				if(sastype==3):
					if(lin[0]=='Electron'):
						trg=lin[3]
						ftrg=locale.atof(trg)
						keep.append(i) 
						rgarray.append(ftrg)
						if(ftrg<rglowcut):
							rgcl.append(1.0)
						else:
							rgcl.append(0.0)
						if(ftrg>rghighcut):
							rgch.append(1.0)
						else:
							rgch.append(0.0)

	if(len(keep)==0):
		message='did not find any profiles that met the input criterion'
		message+=' :  stopping here'
		print_failure(message,txtOutput)
		return
		
	first=0 ; norm=1.0 ; #nh=1  # number of header lines
	for i in xrange(ns):
		if(i in keep):
			if(sastype==1):
				sinf=open(log[i]+'.iq','r').readlines()
				nh=0
			if(sastype==2 or sastype==3):
				sinf=open(log[i]+'.int','r').readlines()
				nh=1
			sln=len(sinf)
			jj=0 
			if(first==0):
				allspec=numpy.zeros((len(keep),(sln-nh)),numpy.float)
				spectra1=numpy.zeros((sln-nh),numpy.float)
				spectra2=numpy.zeros((sln-nh),numpy.float)
				spectra3=numpy.zeros((sln-nh),numpy.float)
				spectra4=numpy.zeros((sln-nh),numpy.float)
				qvar=numpy.zeros((sln-nh),numpy.float)
				first=1
				print '# READING '+str(ns)+' SAS INT FILES'
			for j in xrange(sln):
				slin=string.split(sinf[j])
				if(j>=nh):
					qval=locale.atof(slin[0])
					qvar[jj]=qval
					val1=locale.atof(slin[1])
					val2=locale.atof(slin[2])
					if(sastype!=1):
						val3=locale.atof(slin[3])
						val4=locale.atof(slin[4])
					#if(sastype==3 and jj==0):
					#	norm=val4
					#elif(jj==0):
					if(jj==0):
						norm=val1
					nval1=io*val1/norm
					nval2=io*val2/norm
					#if(sastype==2):
					if(sastype!=1):
						nval3=io*val3/norm
						nval4=io*val4/norm
						allspec[i][jj]=nval1  
					#elif(sastype==3):
					#	nval3=io*val3/norm
					#	nval4=io*val4/norm
					#	allspec[i][jj]=nval4  
					else:
						allspec[i][jj]=nval1	
	
					spectra1[jj]=spectra1[jj]+nval1
					spectra2[jj]=spectra2[jj]+nval2
					if(sastype!=1):
						spectra3[jj]=spectra3[jj]+nval3
						spectra4[jj]=spectra4[jj]+nval4
					jj=jj+1
	
		if(ns>10):
			if((i+1)%(ns/10)==0):
				print '# '+str(i+1)+' FILES READ :\t'+str(int(numpy.ceil((float(i+1)/float(ns))*100.0)))+' percent done'
				fraction_done = (float(i+1)/float(ns))
                        	progress_string='COMPLETED '+str(i+1)+' of '+str(ns)+' : '+str(fraction_done*100.0)+' % done'
                        	print('%s\n' % progress_string)
                        	report_string='STATUS\t'+str(fraction_done)
                        	txtOutput.put(report_string)
	
       	report_string='STATUS\t'+str(1)
       	txtOutput.put(report_string)

	spectra1=spectra1/len(keep)
	spectra2=spectra2/len(keep)
	if(sastype!=1):
		spectra3=spectra3/len(keep)
		spectra4=spectra4/len(keep)
	
	outfile.write('# AVERAGE SPECTRA : file generated on '+ttxt+'\n')
	for k in xrange(sln-nh):
		if(sastype!=1):
			outfile.write('%f\t%f\t%f\t%f\t%f\n' % (qvar[k],spectra1[k],spectra2[k],spectra3[k],spectra4[k]))
		else:
			outfile.write('%f\t%f\t%f\n' % (qvar[k],spectra1[k],spectra2[k]))
	outfile.close()
	
	x2outfile.write('# file generated on '+ttxt+'\n# structure, X2, Rg, filename\n')
	
	numspec=len(keep)
	bq=numpy.zeros(ngpoints,numpy.float)
	best=numpy.zeros(ngpoints,numpy.float)
	worst=numpy.zeros(ngpoints,numpy.float)
	work=numpy.zeros(ngpoints,numpy.float)
	diff=numpy.zeros(ngpoints,numpy.float)
	spec=numpy.zeros((numspec,ngpoints),numpy.float)
	j=0
	for k in xrange(sln-nh):
		print 'k = ',k,' qvar[k] = ',qvar[k]
	
		if(qvar[k]>=low and qvar[k]<=high+epsilon):
			print 'low = ',low,'high = ',high
			print 'j = ',j
			bq[j]=qvar[k]
			#if(sastype==2):
			if(sastype!=1):
				best[j]=spectra1[k]
			#elif(sastype==3):
			#	best[j]=spectra4[k]
			else:
				best[j]=spectra1[k]
			for i in xrange(numspec):
				spec[i][j]=allspec[i][k]
			j+=1
			print 'qvar[k] = ',qvar[k]
			print 'j = ',j

	if(ngpoints!=j): 
		print 'bug'
		print 'ngpoints = ',ngpoints,' j = ',j
		print 'bug'
		print 'low = ',low,'high = ',high
		message='number of grid points '+str(ngpoints)+' BUG!!'
		message+='number of j points '+str(j)+' BUG!!'
		message+=' :  stopping here'
		print_failure(message,txtOutput)

	print '# ANALYZING SAS SPECTRA'	
	search=1 ; nsteps=0
	bestspectra=0
	x2array=[]
	best_tally=[]; worst_tally=[] ; goal_tally=[] ; x2rg_tally=[] ; avg_tally=[] ; avg_sum=[]

	vrarray = []

	if(plotflag == 1):
		graph = Gnuplot.Gnuplot(debug=1)
		graph.clear()
		graph('set title "SAS Spectra"')
		graph('set logscale x')
		graph('set logscale y')
		graph.xlabel('Q (1/Angstrom)')
		graph.ylabel('I(Q)')
		graph2 = Gnuplot.Gnuplot(debug=1)
		graph2.clear()
		graph2('set title "Chi-Square versus Rg"')
		graph2.xlabel('Rg (Angstrom)')
		graph2.ylabel('Chi-Square')

        ttxt=time.ctime()

        st=''.join(['=' for x in xrange(60)])

        txtOutput.put("\n%s \n" %(st))
        txtOutput.put("DATA FROM RUN: %s \n\n" %(ttxt))
	
	for i in xrange(0,nf):
		j=0
		if(sastype!=1): 
			suf='.xiq'
		else:
			suf='.ciq'
		thisfile=open(spectra+'/spec_'+str(i+1)+suf,'w')
		tspec=numpy.zeros(ngpoints,numpy.float)
		for k in xrange(sln-nh):
			if(qvar[k]>=low and qvar[k]<=high+epsilon):
					tspec[j]=allspec[i][j]
					thisfile.write('%f\t%f\t%f\n' % (qvar[j],tspec[j],0.0))
					j=j+1
		thisfile.close()
		diff=tspec-goal

		vr = numpy.sum(abs(diff/((tspec+goal)/2.0)))
		vrarray.append(vr)

		diff2=diff*diff
		if(reduced_x2 == 1):
			X2=numpy.sum(diff2/sigma2)/(ngpoints-1) # reduced X2 (assuming n = 0)
		elif(reduced_x2 == 0):
			X2=numpy.sum(diff2/sigma2) # chi-square distribution
		elif(reduced_x2 == 2):
			X2=numpy.sum(diff2/goal) # Pearson's X2 test-statistic	
		if(X2<x2lowcut):
			x2cl.append(1.0)
		else:
			x2cl.append(0.0)
		if(X2>x2highcut):
			x2ch.append(1.0)
		else:
			x2ch.append(0.0)
		x2rg_tally.append([rgarray[i],X2])
		x2outfile.write('%d\t%f\t%f\t%s\t%f\n' % (i+1,X2,rgarray[i],log[i],vr))	
		x2outfile.flush()
	
		x2array.append(X2)
		if(i==0):
			bX2=X2
			wX2=X2
			print '#INITIAL CHI-SQUARED = ',X2
			worstfilename=log[i]
			best=tspec.copy()
			worst=tspec.copy()
			avg=tspec.copy()
			bestspectra=i
			worstspectra=i
			bestfilename=log[i]
			for j in xrange(ngpoints):
				best_tally.append([gridq[j],tspec[j]])
				worst_tally.append([gridq[j],tspec[j]])
				goal_tally.append([gridq[j],goal[j]])
				#avg_sum.append([gridq[j],tspec[j]])
				avg_sum.append(tspec[j])
		else:	
			for j in xrange(ngpoints):
				avg_sum[j]=avg_sum[j]+tspec[j]
		if(X2<bX2):
			best=tspec.copy()
			bestspectra=i+1
			print 'new best:  this X2 = ',X2
			best_tally=[]
			for j in xrange(ngpoints):
				best_tally.append([gridq[j],tspec[j]])
			bestfilename=log[i]
			bX2=X2
		if(X2>wX2):
			wX2=X2
			print 'new worst:  this X2 = ',X2
			worst=tspec.copy()
			worst_tally=[]
			worstspectra=i+1
			for j in xrange(ngpoints):
				worst_tally.append([gridq[j],tspec[j]])
			worstfilename=log[i]

	for j in xrange(ngpoints):
		avg_tally.append([gridq[j],avg_sum[j]/nf])

	if(plotflag == 1):
		try:
        		graph.plot(Gnuplot.Data(goal_tally,using='1:2 w p pt 4',title='Goal'),Gnuplot.Data(best_tally,using='1:2 w lp pt 5 lw 2',title='Best'),Gnuplot.Data(worst_tally,using='1:2 w lp pt 3 lw 2',title='Worst'),Gnuplot.Data(avg_tally,using='1:2 w lp pt 6 lw 2',title='Average'))
		except:
			message='error trying plot final goal/best/worst/average data'
			message+=' :  stopping here'
			print_failure(message,txtOutput)

	sasplot=open('./'+runst+'sas_spectra_plot.txt','w')
	sasplot.write('# q-value, goal, best, worst, average\n')
	for j in xrange(ngpoints):
		sasplot.write('%f\t%f\t%f\t%f\t%f\n' % (gridq[j],goal_tally[j][1],best_tally[j][1],worst_tally[j][1],avg_tally[j][1]))		

	sasplot.close()
	x2plot=open('./'+runst+'x2_vs_rg_plot.txt','w')
	nnpoints=len(x2rg_tally)
	for j in xrange(nnpoints):
		x2plot.write('%f\t%f\n' % (x2rg_tally[j][0],x2rg_tally[j][1]))
		x2plot.flush()
	x2plot.close()
	x2outfile.close()	
	
	if(plotflag == 1):
		try:
			graph2.plot(Gnuplot.Data(x2rg_tally,using='1:2 w p pt 4',title=''))
		except:
			message='error trying plot x2 versus Rg data'
			message+=' :  stopping here'
			print_failure(message,txtOutput)


	x2lowfile.write('# file generated on '+ttxt+'\n# structure, X2, weight\n')
	x2highfile.write('# file generated on '+ttxt+'\n# structure, X2, weight\n')
	rglowfile.write('# file generated on '+ttxt+'\n# structure, Rg, weight\n')
	rghighfile.write('# file generated on '+ttxt+'\n# structure, Rg, weight\n')

	for i in xrange(0,nf):
		x2lowfile.write('%i\t%f\t%f\n' % (i+1,x2array[i],x2cl[i]))
		x2highfile.write('%i\t%f\t%f\n' % (i+1,x2array[i],x2ch[i]))
		rglowfile.write('%i\t%f\t%f\n' % (i+1,rgarray[i],rgcl[i]))
		rghighfile.write('%i\t%f\t%f\n' % (i+1,rgarray[i],rgch[i]))

	x2lowfile.close()
	x2highfile.close()
	rglowfile.close()
	rghighfile.close()

	print '\n# BEST SPECTRA = ',bestspectra
	print '# BEST CHI-SQUARED = ',bX2
	print '# WORST SPECTRA = ',worstspectra
	print '# WORST CHI-SQUARED = ',wX2
	bestworstoutfile.write('# FILE GENERATED ON '+ttxt+'\n# gridq,best,goal,diff(best-goal),worst,wdiff(worst-goal)\n')
	bestworstoutfile.write('#bestfilename = %s : bX2 = %f, worstfilename = %s : wX2 = %f\n' % (bestfilename,bX2,worstfilename,wX2))
	for k in xrange(ngpoints):
		bestworstoutfile.write('%f\t%f\t%f\t%f\t%f\t%f\n' % (gridq[k],best[k],goal[k],best[k]-goal[k],worst[k],worst[k]-goal[k]))
	bestworstoutfile.close()

	txtOutput.put( "Data stored in directory: %s\n\n" % ('./'+runst))


	txtOutput.put("PROCESSED %i SAS FILES:\n\n" % (nf))
	txtOutput.put(">> The BEST and WORST SAS spectra are in the file named : %s \n" % (bestworstfile))
	txtOutput.put(">> The AVERAGE SAS spectra is in the file named : %s \n" % (averagefile))
	txtOutput.put(">> Chi-square, Rg, and filename are in the file named : %s \n" % (x2file))
	txtOutput.put("\nBEST SINGLE STRUCTURE IS NUMBER %i WTIH X2 = %f : \n\n\t%s\n\n" % (bestspectra,bX2,bestfilename))
	txtOutput.put("\nWORST SINGLE STRUCTURE IS NUMBER %i WTIH X2 = %f: \n\n\t%s\n\n" % (worstspectra,wX2,worstfilename))
	txtOutput.put("\nDENSITY WEIGHT FILES: \n")
	txtOutput.put(">>>> X2 > %f :\t weights file : %s \n" % (x2highcut,'./'+runst+x2highweightfile))
	txtOutput.put(">>>> X2 < %f :\t weights file : %s \n" % (x2lowcut,'./'+runst+x2lowweightfile))
	txtOutput.put(">>>> Rg > %f :\t weights file : %s \n" % (rghighcut,'./'+runst+rghighweightfile))
	txtOutput.put(">>>> Rg < %f :\t weights file : %s \n" % (rglowcut,'./'+runst+rglowweightfile))
	txtOutput.put("\n%s \n" %(st))
        
	if(plotflag == 1):
		wait('\n')
  	else:
  		time.sleep(2)

	print 'CHI_SQUARE IS DONE'

	return

