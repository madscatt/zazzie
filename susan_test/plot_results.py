import os
import glob
import matplotlib.pyplot as plt
import numpy 

count = 0

fig, (ax0,ax1) = plt.subplots(nrows=2)
fig.tight_layout(pad=3.0)
ax0.set_xscale("log")
ax0.set_yscale("log", nonpositive='clip')
ax0.set_title("subunit data", fontsize=11, fontweight="medium")
ax0.set_ylabel('I')
#ax0.set_xlabel('q')
ax1.set_title("inter-subunit data", fontsize=11, fontweight="medium")
ax1.set_ylabel('I')
ax1.set_xlabel('q')

path = './run_0/multi_component_analysis/decomposition/'
for fname in glob.glob(path+'i*run_0.dat'):
    count += 1 
    print('fname: ', fname)
   # data = open(fname,'r').readlines()
    data = numpy.loadtxt(fname)
#    print('data[:,0] = ', data[:,0])

    x = data[:,0]
    y = data[:,1]
    yerror = data[:,2]

    if count == 2:
        ax0.errorbar(x, y, yerr=yerror, fmt='bo', markersize = 2, label = 'I11')
    elif count == 3:
        ax0.errorbar(x, y, yerr=yerror, fmt='ro', markersize = 2, label = 'I22')
    elif count == 1:
        ax1.errorbar(x, y, yerr=yerror, fmt='ko', markersize = 2, label = 'I12') 
# ylim must be set after errorbar to allow errorbar to autoscale limits
#        ax0.set_ylim(bottom=0.00001, top=1.0)
    
'''
#plots
    fig, (ax0,ax1) = plt.subplots(nrows=2)

    ax0.set_xscale("log", nonpositive='clip')
    ax0.set_yscale("log", nonpositive='clip')
    ax0.set(title='Errorbars go negative')
    ax0.errorbar(scattering_data[0], y, xerr=0.1 * x, yerr=5.0 + 0.75 * y)
# ylim must be set after errorbar to allow errorbar to autoscale limits
    ax4.set_ylim(bottom=0.1)

    ax0.errorbar(fraction_d2o, square_root_izero, yerr=square_root_izero_error, fmt='bo', markersize = 2, label = "data")
    ax0.plot(fraction_d2o,ycalc,'r--', label = "calc")
    ax0.set_title('data')
    ax0.set_ylabel('sqrt[I(0)/c]')
    ax0.legend(loc="upper right")    
#
    ax1.plot(fraction_d2o, diff, 'ko', markersize = 5)
    ax1.set_title('residuals')
    ax1.set_xlabel('fraction D2O')
    ax1.set_ylim([-0.5, 0.5])
#    plt.show()

                      

    plt.loglog(x,y)
'''
ax0.legend(loc="upper right") 
ax1.legend(loc="upper right")    
plt.show()

#plt.ylabel('I(q)')
#plt.xlabel('q (1/angstrom)')
#plt.show()
#print 'count = ', count

    #ax.plot(x, y)

    #ax.set_xlabel('X',fontsize=16,fontweight="bold")
    #ax.set_ylabel('Y',fontsize=16,fontweight="bold")
    #plt.savefig(os.getcwd()+fname+'.pdf',figsize=(5,5),dpi=600)
