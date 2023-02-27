import sasmol.sasmol as sasmol
import sys,os, glob, locale, numpy

def get_grid_number(coor,delta):
    idx = int(coor[0]/delta)
    idy = int(coor[1]/delta)
    idz = int(coor[2]/delta)
    return (idx,idy,idz)

def get_number_new_grids(coors,grid_set,delta):
    number_new_grids = 0
    for coor in coors:
        grid_number = get_grid_number(coor,delta)
        if grid_number not in grid_set:
            number_new_grids += 1
            grid_set.add(grid_number)
    return number_new_grids

def converge_real_space(mol, folder):
    grid_set = set([])
    delta = 5.0
    list_number_new_grids = []
    list_number_occupied_grids = []
    number_occupied_grids = 0
    error, mask = mol.get_subset_mask('name[i]=="CA"') 
    file_out = os.path.join(folder,"number_of_occupied_grids_real_space.txt")
    fout = open(file_out,'w')
    fout.write("#frame_number, number_of_occupied_grids\n")
    for nf in xrange(mol.number_of_frames()):
        error, coors = mol.get_coor_using_mask(nf,mask)
        number_new_grids = get_number_new_grids(coors[0],grid_set,delta)
        number_occupied_grids += number_new_grids
        list_number_new_grids.append(number_new_grids)
        list_number_occupied_grids.append(number_occupied_grids)
        fout.write("%d %d\n"%(nf+1, number_occupied_grids))
    fout.close()

def converge_sas_space(runname, folder):
    Iq_all=[]
    steps = []
    for f in glob.glob(os.path.join(folder,runname+'*[0-9].iq')):
        steps.append(locale.atoi(f[-8:-3]))
        Iq=[]
        for line in open(f).readlines()[1:]:
            words=line.split()
            Iq.append(locale.atof(words[1]))
        Iq_all.append(Iq)
    Iq_all = numpy.array(Iq_all)
    Iq_low = numpy.amin(Iq_all,axis=0)
    Iq_high = numpy.amax(Iq_all,axis=0)
    Nspec = Iq_all.shape[0]
    Nq = Iq_all.shape[1]
    Ngrid = 10

    grid = numpy.zeros((Nq,Ngrid+1))
    number_of_occupied_grids = 0
    list_number_of_new_grids = []
    list_number_of_occupied_grids = []
    file_out = os.path.join(folder,"number_of_occupied_grids_sas_space.txt")
    fout = open(file_out,'w')
    fout.write("#sas_profile_number, number_of_occupied_grids\n")
    for nspec in xrange(Nspec):
        number_of_new_grids = 0
        for nq in xrange(Nq):
            if Iq_high[nq]==Iq_low[nq]: continue
            n = int(Ngrid*((Iq_all[nspec,nq]-Iq_low[nq])/(Iq_high[nq]-Iq_low[nq])))
            if grid[nq,n]: continue
            else:
                grid[nq,n] = 1
                number_of_new_grids += 1
        number_of_occupied_grids += number_of_new_grids
        list_number_of_new_grids.append(number_of_new_grids)
        list_number_of_occupied_grids.append(number_of_occupied_grids)
        fout.write("%d %d\n"%(nspec+1, number_of_occupied_grids))
    fout.close()
