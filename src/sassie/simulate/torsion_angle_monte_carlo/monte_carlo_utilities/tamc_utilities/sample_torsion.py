import sys
import math,numpy,random

import sasmol.sasmol as sasmol
import sasmol.sasmath as sasmath

def calc(other_self,itheta,theta,parm,beta,nonbondflag,seed_object):
    '''
    modified from the energy module in simulate package
    ## @NOTE to ZHL: is this ok?
    '''
    log = other_self.log
    pgui = other_self.run_utils.print_gui

    #pgui("initial theta: %f...trial theta: %f"%(itheta,theta))

    search=True; vdi=1E10 ; vdf=1E10
    pi=numpy.pi
    angle_to_radians=pi/180.0
    thetarad=theta*angle_to_radians
    ithetarad=itheta*angle_to_radians
    number_of_parm = len(parm)/3

    iv = 0.0 ; v = 0.0
    for i in xrange(number_of_parm):
        offset = i*3
        k_ang = float(parm[0+offset])
        n_ang = float(parm[1+offset])
        delta_ang = (float(parm[2+offset]))*angle_to_radians

        iv += k_ang*(1.0+math.cos(n_ang*ithetarad-delta_ang))
        v += k_ang*(1.0+math.cos(n_ang*thetarad-delta_ang))

    #pgui("initial energy: %f...trial energy: %f"%(iv,v))

    if(nonbondflag==0):
        if(v<iv):
            #pgui("...lower energy automatically accepted")
            search=False
            vdi=iv
            vdf=v
        else:
            #pgui("...higher energy decide by Metropolis sampling...")
            kcal2kj=4.184 # 1 kcal = 4.184 kJ
            vkjmol=v*kcal2kj
            ivkjmol=iv*kcal2kj
            kjmol2j=1.660540E-21
            vj=vkjmol*kjmol2j
            ivj=ivkjmol*kjmol2j
            earg=(-beta*(vj-ivj))
            if(seed_object != -1):
                r = seed_object.rand()
            else:
                r=random.random()
            arg=numpy.exp(earg)
            if(r<arg):
                #pgui("...accepted!")
                search=False
                vdi=iv
                vdf=v
            else:
                #pgui("...rejected!")
                pass

        return search,vdi,vdf
    else:
        vdi=iv
        vdf=v
        return vdi,vdf

def rotate_dihedral(other_self, group_molecule, group_number, pivot_number, itheta):
    '''
    rotate dihedral
    '''

    result = 1

    mvars = other_self.mvars
    pvars = mvars.pvars

    coor = group_molecule.coor()
    mvars = other_self.mvars
    main_pivots_masks = pvars.main_pivots_masks[group_number]
    main_pivots_post_masks = pvars.main_pivots_post_masks[group_number]
    group_flexible_mask = other_self.group_flexible_masks[group_number]

    frame = 0 ## @NOTE to ZHL: should we pass frame from outside?

    active_pivot_mask = main_pivots_masks[pivot_number]
    #err, active_pivot_coor=group_molecule.get_coor_using_mask(frame, active_pivot_mask)
    tmp_group_flexible_molecule = sasmol.SasMol(0)
    error = group_molecule.copy_molecule_using_mask(tmp_group_flexible_molecule,group_flexible_mask,0)
    err, active_pivot_coor=tmp_group_flexible_molecule.get_coor_using_mask(frame, active_pivot_mask) ## @NOTE to ZHL: temporary solution for a hack
    masks_to_be_rotated = [main_pivots_post_masks[pivot_number]] ## @NOTE to ZHL: need a comprehensive way to build the to-be-rotated masks, including the post-internal and post-external. Hacked currently.
    active_pivot_coor = active_pivot_coor[0] ## @NOTE to ZHL: how do I know this is a forward or backward rotation? Append [::-1] for backward rotation.

    theta=itheta*(math.pi/180.0)

    v=numpy.zeros(3,numpy.float)
    tee=numpy.identity(4,numpy.float)

    r1 = active_pivot_coor[0,:]
    r2 = active_pivot_coor[1,:]
    r3 = active_pivot_coor[2,:]
    r4 = active_pivot_coor[3,:]

    v[0]=r3[0]-r2[0]
    v[1]=r3[1]-r2[1]
    v[2]=r3[2]-r2[2]

    tee[0][3]=r2[0]
    tee[1][3]=r2[1]
    tee[2][3]=r2[2]

    try:
        lv=math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])
        if(lv > 0.0):
            nv=v/lv
        else:
            result = 0
            print 'v = ',v
            print 'lv = ',lv
            return result
    except:
        result = 0
        print 'v = ',v
        print 'lv = ',lv
        return result


    vx=nv[0] ; vy=nv[1] ; vz=nv[2] ; vx2=vx*vx ; vy2=vy*vy ; vz2=vz*vz
    cot=math.cos(theta) ; sit=math.sin(theta)
    r=numpy.zeros((4,4),numpy.float)

    r[0][0]=cot+(1.0-cot)*vx2 ; r[0][1]=(1.0-cot)*vx*vy-sit*vz ; r[0][2]=(1.0-cot)*vx*vz+sit*vy ; r[0][3]=0.0
    r[1][0]=(1.0-cot)*vx*vy+sit*vz ; r[1][1]=cot+(1.0-cot)*vy2 ; r[1][2]=(1.0-cot)*vy*vz-sit*vx ; r[1][3]=0.0
    r[2][0]=(1.0-cot)*vx*vz-sit*vy ; r[2][1]=(1.0-cot)*vy*vz+sit*vx ; r[2][2]=cot+(1.0-cot)*vz2 ; r[2][3]=0.0
    r[3][0]=0.0 ; r[3][1]=0.0 ; r[3][2]=0.0 ; r[3][3]=1.0

    itee = numpy.matrix(tee).I
    ir = numpy.matrix(r)*itee
    br = tee*ir

    for mask in masks_to_be_rotated:
        for i in xrange(len(mask)):
            if not mask[i]:
                continue
            ta=coor[frame][i]
            far=ta.tolist()
            far.append(1.0)
            nta=numpy.array(far)
            nta.shape=(-1,4)
            p=br*numpy.matrix(numpy.transpose(nta)) ; p=numpy.array(p)
            p = numpy.transpose(p)[:,:3]
            coor[frame,i,:] = numpy.transpose(p[0])

    #group_molecule.write_pdb('test.pdb',0,'w'); exit(0) ## @NOTE to ZHL: hardwired output for testing purpose

    return result


def measure(other_self, group_molecule, group_number, pivot_number):

    mvars = other_self.mvars
    pvars = mvars.pvars
    log = other_self.log
    pgui = other_self.run_utils.print_gui

    main_pivots_mask = pvars.main_pivots_masks[group_number][pivot_number]
    coor = group_molecule.coor()

    log.debug('in measure')

    ##pgui('Indices of selected main_pivots_mask = '+str(numpy.nonzero(main_pivots_mask)[0]))

    ind=numpy.nonzero(main_pivots_mask*numpy.arange(1,len(main_pivots_mask)+1))[0]
    this_frame_coor = coor[0,:,:]

    lcoor=numpy.take(this_frame_coor[:,:],ind,0)

    angle=sasmath.dihedral_angle(lcoor[0,:],lcoor[1,:],lcoor[2,:],lcoor[3,:])

    return angle

def check_angle(current_theta,trial_delta_theta):
    '''
    Method checks the trial_delta_theta such that it is between -180 and 180 degrees
    '''

    sum_theta = current_theta + trial_delta_theta

    if (sum_theta > 180.0):
        sum_theta = sum_theta % 180.0
    elif(sum_theta < -180.0):
        sum_theta = abs((sum_theta % 180.0))
    """
    ## @NOTE to ZHL: I don't think the above code is correct, and following is my revision:
    while sum_theta >= 180.0:
        sum_theta -= 360.0
    while sum_theta < -180.0:
        sum_theta += 360.0
    """

    return sum_theta

def find_angle(other_self, group_molecule, group_number, pivot_number):

    mvars = other_self.mvars
    pvars = mvars.pvars
    log = other_self.log
    pgui = other_self.run_utils.print_gui

    main_pivots_parameter = pvars.main_pivots_parameters[group_number][pivot_number]
    kb=1.380658E-23 # J/K
    beta=1.0/(mvars.temperature*kb)
    nonbondflag = mvars.nonbondflag
    seed_object = mvars.seed_object

    log = other_self.log
    log.debug('in find angle')

    current_theta = measure(other_self, group_molecule, group_number, pivot_number)

    search = True

    while search:

        if(seed_object != -1):
            trial_delta_theta = mvars.delta_theta_array[group_number] * (2.0 * seed_object.rand() - 1)
        else:
            trial_delta_theta = mvars.delta_theta_array[group_number] * (2.0 * random.random() - 1)
        new_theta = check_angle(current_theta,trial_delta_theta)

        search,vdi,vdf = calc(other_self,current_theta, new_theta, main_pivots_parameter, beta, nonbondflag, seed_object)

    #pgui("Finally, dihedral energy: %f -> %f"%(vdi,vdf))
    return vdi, vdf, trial_delta_theta

def sample_torsion(other_self, group_number, fout_phi_psi=None):
    '''
    method of TAMC torsion sampling
    '''

    mvars = other_self.mvars
    mcvars = other_self.mcvars
    pvars = mvars.pvars
    log = other_self.log

    group_molecule = other_self.group_molecules[group_number]

    #pgui = other_self.run_utils.print_gui

    log.debug('in sample torsion')

    number_of_pivots = len(pvars.main_pivots_parameters[group_number])
    pivot_number = random.randint(0, number_of_pivots-1) ## @NOTE to ZHL: bug in random_seed_object

    ##pgui('main_pivots_parameters[group_number] = '+str(mvars.main_pivots_parameters[group_number]))
    ##pgui('Selected pivot parameters = '+str(pvars.main_pivots_parameters[group_number][pivot_number]))

    original_phi = measure(other_self, group_molecule, group_number, 1)

    vdi,vdf,trial_delta_theta = find_angle(other_self, group_molecule, group_number, pivot_number)
    rotate_dihedral(other_self, group_molecule, group_number, pivot_number, trial_delta_theta)

    final_phi = measure(other_self, group_molecule, group_number, 1)
    #print ('Phi: %f -> %f\n'%(original_phi,final_phi))

    mcvars.trial_accepted = True

    if fout_phi_psi != None:
        phi = measure(other_self, group_molecule, group_number, 1)
        psi = measure(other_self, group_molecule, group_number, 2)
        fout_phi_psi.write("%f %f\n"%(phi,psi))

    return

