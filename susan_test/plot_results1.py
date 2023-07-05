import os
import glob
import matplotlib.pyplot as plt
import numpy
import read_data_file as read_data_file 

count = 0

fig, (ax0,ax1) = plt.subplots(nrows=2)
fig.tight_layout(pad=3.0)
ax0.set_xscale("log")
ax0.set_yscale("log", nonpositive='clip')
ax0.set_title("rescaled CV data", fontsize=11, fontweight="medium")
ax0.set_ylabel('I')
#ax0.set_xlabel('q')
ax1.set_xscale("log")
ax1.set_yscale("log", nonpositive='clip')
ax1.set_title("rescaled CV data", fontsize=11, fontweight="medium")
ax1.set_ylabel('I')
ax1.set_xlabel('q')

path = './run_0/multi_component_analysis/decomposition/'
odata = ['0.dat','10.dat','20.dat','40.dat','80.dat','90.dat','100.dat']
rdata = ['0_rescaled.dat','10_rescaled.dat','20_rescaled.dat','40_rescaled.dat','80_rescaled.dat','90_rescaled.dat','100_rescaled.dat']
cdata = ['0_calculated.dat','10_calculated.dat','20_calculated.dat','40_calculated.dat','80_calculated.dat','90_calculated.dat','100_calculated.dat']
for fname in rdata:
    count += 1 
    print('fname: ', path+fname)
    data = numpy.loadtxt(path+fname)
#    print('data[:,0] = ', data[:,0])

    x = data[:,0]
    y = data[:,1]
    yerror = data[:,2]

    if count == 1:
        ax0.errorbar(x, y, yerr=yerror, fmt='bo', markersize = 2, label = '0r')
    elif count == 2:
        ax0.errorbar(x, y, yerr=yerror, fmt='ro', markersize = 2, label = '10r')
    elif count == 3:
        ax0.errorbar(x, y, yerr=yerror, fmt='ko', markersize = 2, label = '20r')
    elif count == 4:
        ax0.errorbar(x, y, yerr=yerror, fmt='go', markersize = 2, label = '40r')
    elif count == 5:
        ax1.errorbar(x, y, yerr=yerror, fmt='bo', markersize = 2, label = '80r')
    elif count == 6:
        ax1.errorbar(x, y, yerr=yerror, fmt='ro', markersize = 2, label = '90r')
    elif count == 7:
        ax1.errorbar(x, y, yerr=yerror, fmt='ko', markersize = 2, label = '100r') 
# ylim must be set after errorbar to allow errorbar to autoscale limits
#        ax0.set_ylim(bottom=0.00001, top=1.0)

count = 0
for item in odata:
    count += 1 
    print('fname: ', item)

    x = []
    y = []
    yerr = []
    
    numpts, x, y, yerr = read_data_file.read_file(item)

    if count == 1:
        ax0.errorbar(x, y, yerr=yerror, fmt='bx', markersize = 2, label = '0')
    elif count == 2:
        ax0.errorbar(x, y, yerr=yerror, fmt='rx', markersize = 2, label = '10')
    elif count == 3:
        ax0.errorbar(x, y, yerr=yerror, fmt='kx', markersize = 2, label = '20')
    elif count == 4:
        ax0.errorbar(x, y, yerr=yerror, fmt='gx', markersize = 2, label = '40')
    elif count == 5:
        ax1.errorbar(x, y, yerr=yerror, fmt='bx', markersize = 2, label = '80')
    elif count == 6:
        ax1.errorbar(x, y, yerr=yerror, fmt='rx', markersize = 2, label = '90')
    elif count == 7:
        ax1.errorbar(x, y, yerr=yerror, fmt='kx', markersize = 2, label = '100') 
        ax0.set_ylim(bottom=0.000001, top=1.0)    
        ax1.set_ylim(bottom=0.000001, top=1.0)    
'''
#plots
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
