import sys, string, locale,os
import numpy
import math

import sasmol.sasmol as sasmol

sys.path.append('./')

def get_nh_vector(self):
    '''
   Functionality: Find N-H vector from an input pdb file
                  => Need to include general atom pair given in RDC file header in future
   Input: mol -- sasmol object 
          Note that input pdb should pass pdbrx
   Return: N-H unit vector in residues except N terminal
    '''

    mvars = self.mvars
    mol = sasmol.SasMol(0)
    mol.read_pdb(mvars.pdbfile)
    natoms = mol.natoms()
    coor = mol.coor()
    numres_vnh = mol.number_of_resids()-1
    local_nh_vector_coor = numpy.zeros((numres_vnh,4))

    ''' Check whether atom numbering follows CHARMM convention '''

    count = 0
    second_resid = mol.resids()[1]
 
    for i in xrange(natoms):
        if (mol.resid()[i] == second_resid):
            atom_name_i = mol.name()[i]
            atom_name_j = mol.name()[i+1]
            if (atom_name_i == "N" ):
                try:
                    if(atom_name_j == "HN"):
                           print (">> input pdb seems to pass pdbrx.")
                except: 
                    message = 'Please consider to run pdbrx for correct atom naming'
                    pgui(message)
                    sys.exit(1)

    count = 0 
        
    for i in xrange(natoms-1):
        atom_name_i = mol.name()[i]
        atom_name_j = mol.name()[i+1]
        if (atom_name_i == "N" ):
           if (atom_name_j =="HN"):
               resnumber = mol.resid()[i]
               xx = coor[0,i+1,0] - coor[0,i,0]
               yy = coor[0,i+1,1] - coor[0,i,1]
               zz = coor[0,i+1,2] - coor[0,i,2]
               rr = numpy.sqrt(xx*xx + yy*yy + zz*zz)
               xx = xx/rr
               yy = yy/rr
               zz = zz/rr
           
               local_nh_vector_coor[count][0] = resnumber
               local_nh_vector_coor[count][1] = xx
               local_nh_vector_coor[count][2] = yy
               local_nh_vector_coor[count][3] = zz

               count += 1

    return local_nh_vector_coor
 
def read_in_data(self):

    mvars = self.mvars

    ''' read in RDC data '''
   
    local_rdc = numpy.zeros((len(mvars.rdc),2))
 
    count = 0
    for line in mvars.rdc:
        lin = string.split(line)
        local_rdc[count][0] = locale.atof(lin[0])
        local_rdc[count][1] = locale.atof(lin[1])        
        count += 1

    mvars.rdc = local_rdc
   
    ''' read in vnh vector data '''
   
#    local_nh_vector_coor = numpy.zeros((len(mvars.nh_vector_coor),4))
 
#    count = 0
#    for line in mvars.nh_vector_coor:
#        lin = string.split(line)
#        local_nh_vector_coor[count][0] = locale.atof(lin[0])
#        local_nh_vector_coor[count][1] = locale.atof(lin[1])        
#        local_nh_vector_coor[count][2] = locale.atof(lin[2])        
#        local_nh_vector_coor[count][3] = locale.atof(lin[3])        
#        count += 1

    mvars.nh_vector_coor = get_nh_vector(self)
    n_vnh = len(mvars.nh_vector_coor)

    local_active_residues = []
    count = 0
    for line in mvars.active_residues:
       count += 1
       lin = string.split(line)
       for residue in lin:
           local_active_residues.append(int(locale.atof(residue)))

    mvars.active_residues = local_active_residues

#    print "===== Active residue list =========\n"
#    print mvars.active_residues

    ''' set up random seed '''

    if mvars.mcon:
        if mvars.seed[0] == 1:
            from numpy.random import RandomState
            mvars.seed_object = RandomState(mvars.seed[1])
        else:
            mvars.seed_object = -1

    return

def prep2(self):

    mvars = self.mvars

    ''' no sorting is done here '''
    ''' not checking for NaN '''

    ''' Input variables '''

    di_input = mvars.rdc
    vnh_input = mvars.nh_vector_coor
    reslist_input =mvars.active_residues

    di_len = len(di_input) 
    vnh_len = len(vnh_input)
    reslist_len = len(reslist_input)

    len_max = max(di_len,vnh_len)
    len_max = max(len_max,reslist_len)

    local_di = numpy.zeros((len_max,2))
    local_vnh = numpy.zeros((len_max,4))
    local_reslist = numpy.zeros((len_max))

    local_di_2 = local_di
    local_vnh_2 = local_vnh

    count = 0
    for ii in range(0,vnh_len):
        for jj in range(0,di_len):
            ind = numpy.where(vnh_input[ii][0] == di_input[jj][0])
            if ind[0].size != 0:
               local_vnh[count][0] = di_input[jj][0]
               local_vnh[count][1:] = vnh_input[ii][1:]
               local_di[count][0] = di_input[jj][0]
               local_di[count][1] = di_input[jj][1]
               count = count + 1
    count = 0

    for ii in range(0,len_max):
        for jj in range(0,reslist_len):
            ind = numpy.where(local_di[ii][0] == reslist_input[jj])
            if ind[0].size != 0:
               local_di_2[count][0] = reslist_input[jj]
               local_di_2[count][1] = local_di[ii][1]
               local_vnh_2[count][0] = reslist_input[jj]
               local_vnh_2[count][1:] = local_vnh[ii][1:]
               local_reslist[count] = reslist_input[jj]
               count = count + 1 

    di_output = numpy.zeros((count,2))
    vnh_output = numpy.zeros((count,4))
    reslist_output = numpy.zeros((count))

    for ii in range(0,count):
        di_output[ii][:] = local_di_2[ii][:]
        vnh_output[ii][:] = local_vnh_2[ii][:]
        reslist_output[ii] = local_reslist[ii]

    mvars.di = di_output
    mvars.vcoor = vnh_output
    mvars.rlist = reslist_output

    return

def dir_cos(self):

    mvars = self.mvars

    local_input = mvars.vcoor
    local_input_len = len(local_input)
    local_output = numpy.zeros((len(local_input),3))
 
    for ii in range(0,local_input_len):
         xx = local_input[ii][1]
         yy = local_input[ii][2]
         zz = local_input[ii][3]
         local_output[ii][0] = math.acos(xx)
         local_output[ii][1] = math.acos(yy)
         local_output[ii][2] = math.acos(zz)

    mvars.dircosangle = local_output
  
    return 

def extract_par(m):

    a = m[1][1]
    b = m[2][2]
    c = m[0][1]
    d = m[0][2]
    e = m[1][2]

    output = numpy.array([a,b,c,d,e])
    return output

def gen_a(self):

    mvars = self.mvars

    tol = 0.01
   
    local_input = mvars.dircosangle 
    local_input_len = len(local_input)
    local_output = numpy.zeros((len(local_input),5))

#    local_output = numpy.ones((local_input_len,local_input_len))
    for ii in range(0,local_input_len):
        phix = math.cos(local_input[ii][0])
        phiy = math.cos(local_input[ii][1])
        phiz = math.cos(local_input[ii][2])

        local_output[ii][0] = phiy*phiy - phix*phix
        local_output[ii][1] = phiz*phiz - phix*phix
        local_output[ii][2] = 2.*phix*phiy
        local_output[ii][3] = 2.*phix*phiz
        local_output[ii][4] = 2.*phiy*phiz

    mvars.mat_a = local_output

#    print local_output

    mvars.con_num = numpy.linalg.cond(local_output)

#   test scipy condition number

 
    ''' single value decomposition '''
    ''' numpy svd: u * t * v = a  where t is diagonal element of T matrix'''
    ''' matlab svd: u,t,v = svd(a) while u*t*v_transpose = a  ''' 

    local_u,local_t, local_v = numpy.linalg.svd(local_output,full_matrices=True) 
    local_t_len = len(local_t)


    local_t_mat = numpy.zeros((local_input_len,local_t_len))

    for ii in range(0,local_t_len):
        local_t_mat[ii][ii] = local_t[ii]


    v_trans = numpy.transpose(local_v)
    u_trans = numpy.transpose(local_u)


    mvars.mat_t = local_t_mat
    mvars.mat_v_trans = v_trans

    local_s = numpy.zeros((local_t_len,local_input_len))

    local_max_t = numpy.max(numpy.max(local_t))
    tol_t = tol*local_max_t

    for ii in range(0,5):
        if local_t_mat[ii][ii] > tol_t:
           local_s[ii][ii] = 1./local_t_mat[ii][ii] 
        else:
           local_s[ii][ii] = 0.

    mat_a_inv = numpy.matmul(v_trans,local_s)

    mvars.mat_a_inv = numpy.matmul(mat_a_inv,u_trans)

    return


def x_tensor3(a_inv,d_eff):

    local_a_inv = a_inv
    local_di = d_eff 

    shape_a_inv = local_a_inv.shape
    len_a_inv = shape_a_inv[1]

    local_q_v = numpy.zeros((5))

    local_s = numpy.zeros((3,3))
    for ii in range(0,5):
        for jj in range(0,len_a_inv):
            local_q_v[ii] = local_q_v[ii] + local_a_inv[ii][jj]*local_di[jj][1]

    local_s[0][0] = -local_q_v[0] - local_q_v[1]
    local_s[1][1] = local_q_v[0]
    local_s[2][2] = local_q_v[1]
    local_s[0][1] = local_q_v[2]
    local_s[0][2] = local_q_v[3]
    local_s[1][2] = local_q_v[4]
    local_s[2][0] = local_s[0][2]
    local_s[2][1] = local_s[1][2]
    local_s[1][0] = local_s[0][1]

### Eigen value, eigen vectors

    local_d, local_v = numpy.linalg.eig(local_s)

    local_v = numpy.real(local_v)
    local_d = numpy.real(local_d)
     
    index = numpy.argsort(numpy.abs(local_d))

    local_v = local_v[:,index]

    local_d_copy = local_d

    local_d_copy[:] = local_d[index[:]]

#    mvars.s_tensor = local_s
#    mvars.s_val = local_d_copy
#    mvars.rot = local_v
   
    return local_v, local_d_copy,local_s


def euler_angle(a):

    flag = 0

    d = numpy.array([a,a,a,a,a,a,a,a])

    d[1][2][:] = -d[0][2][:]
    d[2][1][:] = -d[0][1][:]
    d[3][0][:] = -d[0][0][:]
    d[4][1][:] = -d[0][1][:]
    d[4][2][:] = -d[0][2][:]
    d[5][0][:] = -d[0][0][:]
    d[5][2][:] = -d[0][2][:]
    d[6][0][:] = -d[0][0][:]
    d[6][1][:] = -d[0][1][:]
    d[7][:][:] = -d[0][:][:]

    rin = a

    for ii in range(0,8):
        if flag > 0:
           break
        rin = d[ii,:] 
        angle,flag = euler_parm(rin)

    return angle 


def euler2all(a,b,c):

    '''  Assume input vector '''

    z = numpy.ones((4,3))

    zsel = [a,b,c]
    z[0][:] = [a, b, c]
    z[1][:] = [a, b, c+180.]
    z[2][:] = [a+180., 180.-b, 180.-c]
    z[3][:] = [a+180., 180.-b, 360.-c]

    ''' Convert angles to normal range '''

    for ii in range(0,4):
        for jj in range(0,3):
            if z[ii][jj] >= 270.:
               z[ii][jj] = z[ii][jj] -360.
            if z[ii][jj] <= -270.:
               z[ii][jj] = z[ii][jj] +360.
    '''  Select  '''''

    for ii in range(0,4):
        if (z[ii][0] >= 0.) & (z[ii][0] <= 180.) & (z[ii][1] >= 0.) & (z[ii][1] <= 180.) & (z[ii][2] >= 0.) & (z[ii][2] <= 180.): 
            zsel[:] = z[ii][:]                  
      
    return z, zsel

def euler_parm(rin):

#   rin, r: 3x3 matrix

    r = rin

    tol = 1.0e-12
    tol_sin = 1.0e-12

    angle = numpy.zeros((3))
    zflag = 0
    negative_flag = 1

    if negative_flag > 0:
       if numpy.linalg.det(r) < 0:
          r[2][:] = -rin[2][:]

    t = math.acos(r[2][2])

    sin_t = math.sin(t)

    if abs(sin_t) > tol_sin:
       f1 = math.asin(r[2][1]/sin_t)
       f2 = math.pi - f1
       k1 = math.asin(r[1][2]/sin_t)
       k2 = math.pi - k1

       flag = 0

       flag = choose_euler(r,flag,f1,t,k1,0)
       flag = choose_euler(r,flag,f1,t,k2,1)
       flag = choose_euler(r,flag,f2,t,k1,2)
       flag = choose_euler(r,flag,f2,t,k2,3)

       result = numpy.array([[0,f1,t,k1],[1,f1,t,k2],[2,f2,t,k1],[3,f2,t,k2]])

# This may be redundant

#       for ii in range (0,4):
#           index = numpy.where(flag == result[ii][0])
#

#       print "In Euler_parm index = ", index

#       if index[0].size != 0:
#          idd = index[0][0]
       for jj in range (0,3):
           angle[jj] = 180.*(result[flag][jj+1])/math.pi

       for ii in range(0,3):
          if angle[ii] < 0.:
             angle[ii] = 360. - numpy.abs(angle[ii])
    else:

        if math.cos(t)+1 < tol:
           angle[1] = 180.
           angle[0] = (math.acos(r[0][0]))*180./math.pi
           angle[2] = 0.
        if math.cos(t)-1 < tol:
           angle[1] = 0.
           angle[0] = (math.acos(r[0][0]))*180./math.pi
           angle[2] = 0.
    if (angle[0] != 0.)&(angle[1] != 0.)& (angle[2] != 0.):
       zflag = 1
 
    return angle,flag 

def choose_euler(r,flag,a,b,c,num):

    tol = 1.0e-12

    if flag == 0:
       if (math.cos(a)*math.cos(b)*math.cos(c)-math.sin(a)*math.sin(c)) - r[0][0] < tol:
          if (math.sin(a)*math.cos(b)*math.cos(c)+math.cos(a)*math.sin(c)) - r[0][1] < tol:
             if (-math.cos(a)*math.cos(b)*math.sin(c)-math.sin(a)*math.cos(c)) - r[1][0] < tol:
                if (-math.sin(a)*math.cos(b)*math.sin(c)+math.cos(a)*math.cos(c)) - r[1][1] < tol:
                   if (math.cos(a)*math.sin(b)) - r[2][0] < tol:
                      flag = num  
    
    return flag


def recalc_d(a,b):

    new = numpy.matmul(b,a)

    return new


def calc_r(exp,calc):

    dobs_dcalc = exp - calc
    num = numpy.matmul(dobs_dcalc,dobs_dcalc)
    num = numpy.mean(num)

    denom = 2.*numpy.mean(numpy.matmul(exp,exp))

    r_coeff = math.sqrt(num/denom)

    return r_coeff


def calc_ani(s):

    aniso = (s[0]-s[1])/s[2]

    return aniso


def mc_table(self,a):

    mvars = self.mvars

    num = len(a[:,1])

    output = numpy.zeros((num,2)) 

#   Generate random errors for RDCs 

    ran_num = mvars.seed_object.randn(num,2)

    output[:,1] = a[:,1] + ran_num[:,1] 

    return output


def altens_core(self,app):

    mvars = self.mvars
    self.log.debug('in altens_core')
    pgui = self.run_utils.print_gui

    #mvars.runname = variables['runname'][0]


    mvars.rdc = open(mvars.rdc_input_file,'r').readlines()        
    mvars.active_residues = open(mvars.residue_list_file,'r').readlines() 
 
    read_in_data(self)
 
    prep2(self)

    mvars.max_di = numpy.amax(mvars.di[:,1]) 
    mvars.min_di = numpy.amin(mvars.di[:,1])


    pgui("calculate direction cosine")

    dir_cos(self)

    pgui("Generate matrix A")

    gen_a(self)

    pgui("Single value decomposition")

    eigenvect, s, s_tensor = x_tensor3(mvars.mat_a_inv,mvars.di)

    eigenvect_transpose = numpy.transpose(eigenvect)
    e_ang = euler_angle(eigenvect_transpose)

    all_eulers, select_eulers = euler2all(e_ang[0],e_ang[1],e_ang[2])

    e_ang = select_eulers
    prin_angle = e_ang
    prin_param = s

    ''' extract Sij and recalculate Dij '''

    s_elements = extract_par(s_tensor)

    new_d = recalc_d(s_elements,mvars.mat_a)

#    ds = numpy.array([mvars.rlist,mvars.di,new_d,mvars.di-new_d])

    ''' calculate R factor '''

    corr = numpy.corrcoef(mvars.di[:,1],new_d)

    r = calc_r(mvars.di[:,1],new_d)

    an_param = calc_ani(prin_param)

    ''' write summary '''

    frame = 0
    output_dir = os.path.join(mvars.runname, app)
    output_rdc = mvars.runname+'_calc_rdc_'+str(frame+1).zfill(5)+'.txt'

    fout = open(os.path.join(output_dir, mvars.runname+'_'+str(frame+1).zfill(5))+'.txt', 'w')
    frdc = open(os.path.join(output_dir, mvars.runname+'_calc_rdc_'+str(frame+1).zfill(5))+'.txt','w')

    fout.write("#Altens Information:\n")
    fout.write("1. PDB input file = %s\n"%mvars.pdbfile)
    fout.write("2. RDC input file = %s\n"%mvars.rdc_input_file)
    fout.write("3. residue list file = %s\n"%mvars.residue_list_file)
    fout.write("4. Monte Carlo analysis = %s\n"%mvars.mcon)
    fout.write("5. Calculated RDC file = %s\n"%output_rdc)
    fout.write("6. Correlation ftn of RDC(calc) and RDC(exp) = %.5f\n"%corr[0][1])
    fout.write("7. R factor of calculated RDC = %.5f\n"%r)
    fout.write("8. Alignment tensor \n")
    fout.write(" %.5f %.5f  %.5f\n"%(s_tensor[0][0], s_tensor[0][1], s_tensor[0][2]))
    fout.write(" %.5f %.5f  %.5f\n"%(s_tensor[1][0], s_tensor[1][1], s_tensor[1][2]))
    fout.write(" %.5f %.5f  %.5f\n"%(s_tensor[2][0], s_tensor[2][1], s_tensor[2][2]))
    fout.write("9. Eigenvalues: Sxx, Syy, Szz = %.5f, %.5f, %.5f\n"%(s[0],s[1],s[2]))
    fout.write("10. Euler angles: alpha, beta, gamma = %.5f, %.5f, %.5f\n"%(e_ang[0],e_ang[1],e_ang[2])) 

    fout.close()

    if mvars.mcon:
        pgui("Monte-Carlo analysis")
        mc_test(self,app,prin_param,prin_angle)

    ''' save output rdc files '''
    len_out_rdc = len(new_d)
    frdc.write("Resid          RDC_exp        RDC_calc     delta(Exp-Calc)\n")
    for ii in xrange(len_out_rdc):
        frdc.write("%12.5f %12.5f %12.5f %12.5f \n"%(mvars.di[ii][0], mvars.di[ii][1], new_d[ii],mvars.di[ii][1]-new_d[ii]))
#        print mvars.rlist[ii], mvars.di[ii][1], new_d[ii],mvars.di[ii][1]-new_d[ii]
    frdc.close()

    pgui("Work with figures") 
    ''' 
    Plot the figures below by using output files ...

    Figure 1. mvars.rdc_input_file (bar graph)
              #resid  RDC_exp
              x = column_0 (label: Residue number)
              y = column_1 (label: Experimental RDC)
    Figure 2. frdc (run_0_calc_rdc_00001.txt) 
              #resid  RDC_exp RDC_calc Delta
              x = column_0 (label: Residue number)
              y = column_3 (label: RDCexp - RDCcalc, unit=Hz) 
    Figure 3. frdc (run_0_calc_rdc_00001.txt)
              #resid  RDC_exp RDC_calc Delta
              x = column_1 (label: Experimental RDC) 
              y = column_2 (label: Calculated RDC) 
              Add linear regression
    if mcon: 
    Figure 4. Histogram of simulated alpha,beta,gamma, sxx, syy, szz
              use run_0_mc_trajectory_00001.txt
              #step Sxx Syy Szz alpha beta gamma 
    '''
               
    return


def mc_test(self,app,prin_param,prin_angle):

    ''' Monte Carlo error test '''
    mvars = self.mvars
    nsteps = mvars.number_of_monte_carlo_steps 
    e_set = numpy.ones((nsteps,3))
    s_set = numpy.ones((nsteps,3))
    n_gen = 0
    k = 0

    ''' Consider multi-frame work'''
    frame = 0

    ''' setup trajectory file'''
    output_dir = os.path.join(mvars.runname, app)
    ftr = open(os.path.join(output_dir, mvars.runname+'_'+'mc_trajectory_'+str(frame+1).zfill(5))+'.txt', 'w')
    fout = open(os.path.join(output_dir, mvars.runname+'_'+'mc_analysis_'+str(frame+1).zfill(5))+'.txt', 'w')

    for ii in range(0,nsteps):
        n_gen = n_gen + 1
        d_old = mvars.di
        d_new_mc = mc_table(self,d_old)
        eigenvect_mc, s_mc, s_tensor_mc = x_tensor3(mvars.mat_a_inv,d_new_mc)

        eigenvect_mc_transpose = numpy.transpose(eigenvect_mc)
        e_ang = euler_angle(eigenvect_mc_transpose)
        all_eulers, select_eulers = euler2all(e_ang[0],e_ang[1],e_ang[2])
        e_ang = select_eulers

        if ii == 0:
           e_ang_full = numpy.array([e_ang])
        else:
           e_ang_full  = numpy.append(e_ang_full,numpy.array([e_ang]),axis = 0)

#    ''' check if S principal parameters are meaningful '''

        if ((s_mc[0]<=(prin_param[0]+5.0))&(s_mc[0]>=(prin_param[0]-5.0))&(s_mc[1]<=(prin_param[1]+5.0))&(s_mc[1]>=(prin_param[1]-5.0))& \
           (s_mc[2]<=(prin_param[2]+5.0))&(s_mc[2]>=(prin_param[2]-5.0))&(e_ang[0]<=(prin_angle[0]+60.0))&(e_ang[0]>=(prin_angle[0]-60.0)) \
           &(e_ang[1]<=(prin_angle[1]+60.0))&(e_ang[1]>=(prin_angle[1]-60.0))&(e_ang[2]<=(prin_angle[2]+60.0))&(e_ang[2]>=(prin_angle[2]-60.0))):
           k=k+1

           if ii == 0:
              tab_param = numpy.array([[s_mc[0],s_mc[1],s_mc[2]]]) 
              tab_ang  = numpy.array([[e_ang[0],e_ang[1],e_ang[2]]])
           else:
              tab_param = numpy.append(tab_param,numpy.array([[s_mc[0],s_mc[1],s_mc[2]]]),axis=0)
              tab_ang = numpy.append(tab_ang,numpy.array([[e_ang[0],e_ang[1],e_ang[2]]]),axis=0)
        ftr.write("%12.5f  %12.5f %12.5f  %12.5f  %12.5f  %12.5f\n" \
             %(s_mc[0],s_mc[1],s_mc[2],e_ang[0],e_ang[1],e_ang[2]))
#       test_progress(k,nsteps)

     
    sxx = numpy.mean(tab_param,axis=0)[0]
    std_sxx = numpy.std(tab_param,axis=0)[0]

    syy = numpy.mean(tab_param,axis=0)[1]
    std_syy = numpy.std(tab_param,axis=0)[1]

    szz = numpy.mean(tab_param,axis=0)[2]
    std_szz = numpy.std(tab_param,axis=0)[2]

    res_alpha = numpy.mean(tab_ang,axis=0)[0]
    sd_res_alpha = numpy.std(tab_ang,axis=0)[0]

    res_beta = numpy.mean(tab_ang,axis=0)[1]
    sd_res_beta = numpy.std(tab_ang,axis=0)[1]

    res_gamma = numpy.mean(tab_ang,axis=0)[2]
    sd_res_gamma = numpy.std(tab_ang,axis=0)[2]

    param_mc = [sxx,syy,szz]
    an_param_mc = calc_ani(param_mc)
    error_ani_mc = an_param_mc*math.sqrt(((std_sxx*std_sxx + std_syy*std_syy)/(sxx-syy)/(sxx-syy)) + \
                       (std_szz*std_szz/szz/szz))

    fout.write("========== Monte Carlo error Report ============\n")

    fout.write("1. Number of generated parameters = %s\n"%n_gen)
    fout.write("2. Sxx = %.5f +- %.5f\n"%(sxx,std_sxx))
    fout.write("3. Syy = %.5f +- %.5f\n"%(syy,std_syy))
    fout.write("4. Szz = %.5f +- %.5f\n"%(szz,std_szz))

    fout.write("5. alpha = %.5f +- %.5f\n"%(res_alpha,sd_res_alpha))
    fout.write("6. beta = %.5f +- %.5f\n"%(res_beta,sd_res_beta))
    fout.write("7. gamma = %.5f +- %.5f\n"%(res_gamma,sd_res_gamma))

    fout.write("8. asymmetry parameter = %.5f +- %.5f\n"%(an_param_mc,error_ani_mc))
    
    fout.close()

    return
