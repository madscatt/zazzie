import sys, string, locale
import numpy

sys.path.append('./')

def read_in_data(self):
    mvars = self.mvars

    ''' NOTE: this assumes active residues is a single line of numbers separated by spaces '''

    local_active_residues = []

    #st = locale.atof(string.split(lin[1])[0])

    count = 0
    for line in mvars.active_residues:
        count += 1
        lin = string.split(line)
        for residue in lin:
            local_active_residues.append(int(locale.atof(residue)))

#    print 'count = ',count

#    print 'len local_active_residues = ',len(local_active_residues), '\nlocal_active_residues = ', local_active_residues

    mvars.active_residues = local_active_residues

    ''' read in RDC data '''
   
    local_rdc = numpy.zeros((len(mvars.rdc),2))
 
    count = 0
    for line in mvars.rdc:
        lin = string.split(line)
        local_rdc[count][0] = locale.atof(lin[0])
        local_rdc[count][1] = locale.atof(lin[1])        
        count += 1

    mvars.rdc = local_rdc
#    print 'len rdc = ',len(mvars.rdc), '\nrdc = ', mvars.rdc

    ''' read in NH vector data '''
   
    local_nh_vector_coor = numpy.zeros((len(mvars.nh_vector_coor),4))
 
    count = 0
    for line in mvars.nh_vector_coor:
        lin = string.split(line)
        local_nh_vector_coor[count][0] = locale.atof(lin[0])
        local_nh_vector_coor[count][1] = locale.atof(lin[1])        
        local_nh_vector_coor[count][2] = locale.atof(lin[2])        
        local_nh_vector_coor[count][3] = locale.atof(lin[3])        
        count += 1

#   print 'count = ',count

    mvars.nh_vector_coor = local_nh_vector_coor
#    print 'len nh_vector_coor = ',len(mvars.nh_vector_coor), '\nnh_vector_coor = ', mvars.nh_vector_coor

#    print 'nvc[-1] = ',mvars.nh_vector_coor[-1]
#    print 'nvc[-1][0] = ',mvars.nh_vector_coor[-1][0]

    return

def prep2(self):

    mvars = self.mvars

    ''' no sorting is done here '''
    ''' not checking for NaN '''

    NH = mvars.nh_vector_coor
    nNH = numpy.array(len(mvars.nh_vector_coor))

    nres0 = len(mvars.active_residues)

    ''' prepare and normalize nh_vectors '''

    vNH = numpy.ones((nres0,4))
    vNH[:,0] = mvars.active_residues

    print 'nres0 = ',nres0

## here

    for ii in xrange(0,nres0):
        ind = numpy.where(NH[:][0]==vNH[ii][0])
        if ind[0].size != 0:
            print 'this_vnh = ',vNH[ii][0],'\nind = ',ind,'\nind[0][0] = ',ind[0][0]
            vNH[ii][1:4] = NH[ind[0][0]][1:4]/numpy.sqrt(numpy.sum(NH[ind[0][0]][1:4]*NH[ind[0][0]][1:4]))
            
            print 'in loop: NH[ind[0][0][1:4] = ', NH[ind[0][0]][1:4]
            print 'in loop: vNH[ii] = ', vNH[ii]
            print 'in loop (norm) = ', NH[ind[0][0]][1:4]/numpy.sqrt(NH[ind[0][0]][1:4]*NH[ind[0][0]][1:4])
            print 'in loop (value) = ', NH[ind[0][0]][1:4]
             
    print 'vNH[:][0] = ', vNH[:][0]
    print 'NH[:][0] = ', NH[:][0]
#    print 'vNH = ', vNH
#    print 'NH = ', NH

    return



def altens_core(self, **kwargs):

    mvars = self.mvars
    self.log.debug('in altens_core')
    pgui = self.run_utils.print_gui

    #mvars.runname = variables['runname'][0]

    if kwargs['use_input_files']:

        mvars.rdc = open(mvars.rdc_input_file,'r').readlines()        
        mvars.nh_vector_coor = open(mvars.nh_vector_coordinate_file,'r').readlines()        
        mvars.active_residues = open(mvars.active_residue_list_file,'r').readlines()        
        read_in_data(self)
 
    #mvars.number_of_monte_carlo_steps = variables['number_of_monte_carlo_steps']

    prep2(self)







    return
