import os
import sys
import locale
import string
import sassie.interface.input_filter as input_filter
import sassie.sasmol.sasmol as sasmol
import sassie.simulate.constraints.constraints as constraints
import sassie.simulate.monte_carlo.monomer.overlap as overlap


#used for SASSIE_gui?  Not tested.
def check_psegvariables_exist(psegvariables, glbls):
    error = []
    if psegvariables not in glbls:
        error.append(
            'Please enter the dihedral variables by clicking on the button "click to enter dihedral variables"!')
    return error


def check_segvariables(pdbfile, psegvariables, flpsegname, seglow, seghigh):

    error = []

    allsith = []
    allsnumranges = []
    allsrlow = []
    allsrnum = []
    allmoltype = []

    for i in range(len(psegvariables)):
        allsnumranges.append(psegvariables[i][0])
        allsith.append(psegvariables[i][1])
        allsrlow.append(psegvariables[i][2])
        allsrnum.append(psegvariables[i][3])
        allmoltype.append(psegvariables[i][4])

    '''
	for i in range(len(psegvariables)):
		moltype = allmoltype[i]
		if moltype.strip() not in ['protein','rna']:
			error.append('The molecule type "'+moltype+'" for flexible segment number '+str(i)+' in the flexible segment input fields should be either protein or rna!')
			return error
	'''

    aith = []
    anumranges = []
    arlow = []
    arnum = []
    for i in range(len(allsith)):
        linith = string.split(allsith[i], ',')
        locith = []
        for i in range(len(linith)):
            tith = linith[i]
            try:
                fith = locale.atof(tith)
            except:
                error.append('The angle value "' + tith +
                             '" should be a float type!')
                return error
            if(fith > 180.0 or fith <= 0.0):
                error.append('The angle value "' + tith +
                             '" should be in the range of (0.0,180.0)!')
                return error
            '''
			if(fith>180.0):
				fith=180.0
			elif(fith<0.0):
				fith=0.0
			'''
            locith.append(fith)
        aith.append(locith)

    for i in range(len(allsnumranges)):
        try:
            nr = locale.atoi(allsnumranges[i])
        except:
            error.append('The number of ranges "' + allsnumranges[i] + '" for flexible segment number ' + str(
                i+1) + ' in the flexible segment input fields should be an integer type!')
            return error
        if nr < 1:
            error.append('The number of ranges "' + allsnumranges[i] + '" for flexible segment number ' + str(
                i+1) + ' in the flexible segment input fields should be equal/greater than 1!')
            return error
        anumranges.append(nr)   
#    print 'len(psegvariables), anumranges: ',len(psegvariables), anumranges

    for i in range(len(allsrlow)):
        linrlow = string.split(allsrlow[i], ',')
        linrnum = string.split(allsrnum[i], ',')
        rlow = []
        rnum = []
        for k in range(len(linrlow)):
            try:
                trlow = locale.atoi(linrlow[k])
            except:
                error.append('The low resid "' + allsrlow[i] + '" for flexible segment number ' + str(
                    i+1) + ' in the flexible segment input fields should be an integer array!')
                return error
            rlow.append(trlow)
        for k in range(len(linrnum)):
            try:
                trnum = locale.atoi(linrnum[k])
            except:
                error.append('The number of contiguous residues "' + allsrnum[i] + '" for flexible segment number ' + str(
                    i+1) + ' in the flexible segment input fields should be an integer array!')
                return error
            rnum.append(trnum)
        arlow.append(rlow)
        arnum.append(rnum)

    
    for i in range(len(psegvariables)):
        numranges = anumranges[i]

        #NOT TESTED; test #2 in try statement above will be triggered before getting to this test.
        if numranges < 1:                                                                     
            error.append('number of ranges ' + str(numranges) +
                         ' should be equal or greater than 1!')
            return error

        locvariables = ['resid', 'moltype']
        segname = (flpsegname.split(','))[i].strip()
        value, result = input_filter.get_pdb_complex_stats(
            pdbfile, segname, locvariables)
        if value == 1:
            resid = map(int, result[0])
            # moltype=map(str,result[1])
            moltype = result[1]
#value==0 doesn't seem to occur if pdb file passes tests for valid pdb file
#Not tested, but kept just in case
        elif value == 0:
            error.append(
                'cannot get the pdb statistics for segment ' + segname)
            return error
        number_aa = resid[-1] - resid[0] + 1

        seg_moltype = allmoltype[i]
        if seg_moltype.strip() not in moltype:
            error.append('The molecule type "' + seg_moltype + '" provided for flexible segment number ' +
                         str(i+1) + ' in the flexible segment input fields does not match the pdb file!')
            return error

        lowres1 = seglow[i]
        highres1 = seghigh[i]
        if(lowres1 not in resid):
            error.append('input pdb file, ' + str(pdbfile) + ' does not have low alignment amino acid residue, ' +
                         str(lowres1) + ', range = ' + str(resid[0]) + ' : ' + str(resid[-1]))
            return error
        elif(highres1 not in resid):
            error.append('input pdb file, ' + str(pdbfile) + ' does not have high alignment amino acid residue, ' +
                         str(highres1) + ', range = ' + str(resid[0]) + ' : ' + str(resid[-1]))
            return error
        elif(highres1 - lowres1 < 3):
            error.append(
                'alignment basis is too small (less than 3 points) or low residue > high residue')
            return error

        dtheta = aith[i]
        reslow = arlow[i]
        numcont = arnum[i]
        if(len(dtheta) != numranges):
            error.append('the number of dtheta values does not match the number of ranges, dtheta = ' +
                         str(dtheta) + ' numranges = ' + str(numranges))
            return error
        elif(len(reslow) != numranges):
            error.append('the number of low residue values does not match the number of ranges, lowres = ' +
                         str(reslow) + ' numranges = ' + str(numranges))
            return error
        elif(len(numcont) != numranges):
            error.append('the number of contiguous residues does not match the number of ranges, contiguous = ' +
                         str(numcont) + ' numranges = ' + str(numranges))
            return error

        ''' #Checked above already
		for th in dtheta:
			if(th > 180.0):
				dtheta[th]=180.0
			elif(th < 0.0):
				dtheta[th]=0.0
		'''

        for j in xrange(numranges):
            if (reslow[j] not in resid):
                error.append('Input pdb file, ' + str(pdbfile) + ' does not have low residue amino acid, "' + str(
                    reslow[j]) + '" for segment number ' + str(j+1) + ', range = ' + str(resid[0]) + ' : ' + str(resid[-1]))
                return error
            elif (reslow[j] + numcont[j] not in resid):
                error.append('Input pdb file, ' + str(pdbfile) + ' does not have low+contiguous residue amino acid, "' + str(
                    reslow[j] + numcont[j]) + '" for segment number ' + str(j+1) + ', range = ' + str(resid[0]) + ' : ' + str(resid[-1]))
                return error
            elif (numcont[j] <= 0):
                error.append('The number of contiguous residues "' +
                             str(numcont[j]) + '" should be greater than 0!')
                return error
            elif (not ((lowres1 > (reslow[j] + numcont[j])) or (highres1 < reslow[j]))):
                error.append(
                    'alignment and flexible ranges should not overlap!')
                return error

#if numranges > 1:
        for j in xrange(numranges - 1):
            print 'j, reslow[j]: ', j, reslow[j]
            if(reslow[j] < 2):
                error.append(
                    'low residue can not include the n-terminus, reslow = ' + str(reslow))
                return error
            elif(reslow[j] > reslow[j + 1]):
                error.append(
                    'low residue values must increase from low to high, reslow = ' + str(reslow[j]))
                return error
            elif(reslow[j] + numcont[j] > reslow[j + 1]):
                error.append('low residue values plus number contiguous overlap, reslow = ' +
                             str(reslow[j]) + ' numcont = ' + str(numcont[j]))
                return error
            elif(reslow[j] + numcont[j] > number_aa - 1):       
                error.append('your low residue plus number contiguous exceeds the number of amino acids-1 (' +
                             str(number_aa-1) + '), reslow = ' + str(reslow[j]) + ' numcont = ' + str(numcont[j]))
                return error
        if(reslow[-1] + numcont[-1] > number_aa - 1):             
            error.append('your low residue plus number contiguous exceeds the number of amino acids-1 (' +
                         str(number_aa-1) + '), reslow = ' + str(reslow[-1]) + ' numcont = ' + str(numcont[-1]))
            return error

    return error


def check_complex(variables, psegvariables, **kwargs):

    runname = variables['runname'][0]
    dcdfile = variables['dcdfile'][0]
    path = variables['path'][0]
    pdbfile = variables['pdbfile'][0]
    trials = variables['trials'][0]
    goback = variables['goback'][0]
    temp = variables['temp'][0]

    nsegments = variables['nsegments'][0]
    segbasis = variables['segbasis'][0]
    npsegments = variables['npsegments'][0]
    flpsegname = variables['flpsegname'][0]
    seglow = variables['seglow'][0]
    seghigh = variables['seghigh'][0]

    #basis		= variables['basis'][0]
    #cutoff	    = variables['cutoff'][0]
    lowrg = variables['lowrg'][0]
    highrg = variables['highrg'][0]
    directedmc = variables['directedmc'][0]
    zflag = variables['zflag'][0]
    zcutoff = variables['zcutoff'][0]
    cflag = variables['cflag'][0]
    confile = variables['confile'][0]

    error = []
    error = input_filter.check_name(runname)
    if(error != []):
        return error

    if 'no_file_check' not in kwargs:
        ev, rv, wv = input_filter.check_permissions(path)
        if(ev == 0 or rv == 0 or wv == 0):
            error.append('permission error in input file path ' +
                         path + '  [code = ' + str(ev) + str(rv) + str(wv) + ']')
            if(ev == 0):
                error.append('path does not exist')
            elif(rv == 0):
                error.append('read permission not allowed')
            elif(wv == 0):
                error.append('write permission not allowed')
            return error

        pdbfile = path + '/' + pdbfile

    ev, value = input_filter.check_pdb_dcd(pdbfile, 'pdb')

    if(ev == 0):
        error.append('input pdb file, ' + pdbfile + ', does not exist')
        return error
    if(value == 0):
        error.append('input pdb file, ' + pdbfile +
                     ', is not a valid pdb file')
        return error

    if(trials < 1):
        error.append('trials = ' + str(trials) + '?')
        return error
    elif(goback < 1):
        error.append('goback = ' + str(goback) + '?')
        return error
    elif(temp < 0):
        error.append('use a positive temperature, temperature = ' + str(temp))
        return error
    # elif(cutoff < 1.0):
        #error.append( 'use a larger cutoff value, cutoff = '+str(cutoff))
    #	return error
    elif(zflag != 0 and zflag != 1):
        error.append(
            'zflag == 0 for "no" and 1 for "yes", zflag = ' + str(zflag))
        return error
   # elif(zflag == 1):
   #     if zcutoff < 1.0:
   #         error.append( 'use a larger zcutoff value, zcutoff = '+str(zcutoff))
   #         return error
    elif(cflag != 0 and cflag != 1):
        error.append(
            'cflag == 0 for "no" and 1 for "yes", cflag = ' + str(cflag))
        return error
    elif(cflag == 1):
        locerror = input_filter.check_file_exists(os.path.join(path, confile))
        if len(locerror):
            error.append(locerror)
            return error
        filter_flag = 1
        m1 = sasmol.SasMol(0)
        m1.read_pdb(pdbfile)
        err = constraints.read_constraints(m1, confile, filter_flag)
        if(err != []):
            error.append(err[0])
            return error

    elif(lowrg > highrg):
        error.append('low Rg cutoff is larger than high Rg cutoff, lowrg = ' +
                     str(lowrg) + ' highrg = ' + str(highrg))
        return error
    elif(lowrg < 0 or highrg < 0):
        error.append('Rg cutoffs need to be >= zero, lowrg = ' +
                     str(lowrg) + ' highrg = ' + str(highrg))
        return error
    elif(directedmc < 0):
        error.append(
            'directed Monte Carlo needs to be 0 or a float > 0 (the "goal Rg") ... you entered: ' + str(directedmc))
        return error
    elif(nsegments < 1):
        error.append('the number of segments ' + str(nsegments) +
                     ' should be equal/greater than 1!')
        return error
    elif(npsegments < 1):
        error.append('the number of flexible segments ' +
                     str(npsegments) + ' should be equal/greater than 1!')
        return error
    elif (npsegments > nsegments):
        error.append('the number of flexible segments ' + str(npsegments) +
                     ' should be equal/less than the number of total segments ' + str(nsegments) + '!')
        return error

    segbasis = segbasis.strip()

    if(segbasis.lower() != 'all' and segbasis.lower() != 'heavy' and segbasis.lower() != 'backbone'):
        if(len(segbasis.split(',')) != nsegments):
            error.append('the number of segment basis does not match the number of segments: number of segbasis = ' +
                         str(len(segbasis.split(','))) + ' number of segments = ' + str(nsegments))
            error.append(
                'segment overlap basis entries can be "heavy", "backbone", "all", or a comma delimited list of atom names ... one for each segment')
            return error

    if(len(flpsegname.split(',')) != npsegments):
        error.append('the number of flexible segment names does not match the number of flexible segments: number of flexible segment names = ' +
                      str(len(flpsegname.split(','))) + ' number of flexible segments = ' + str(npsegments))
        return error
    elif(len(seglow) != npsegments):
        error.append('the number of alignment low residues does not match the number of flexible segments: number of alignment low residues = ' +
                     str(len(seglow)) + ' number of  flexible segments = ' + str(npsegments))
        return error
    elif(len(seghigh) != npsegments):
        error.append('the number of alignment high residues does not match the number of flexible segments: number of alignment high residues = ' +
                     str(len(seghigh)) + ' number of  flexible segments = ' + str(npsegments))
        return error

    m1 = sasmol.SasMol(0)
    m1.read_pdb(pdbfile, fastread=True)
    name = m1.name()

    check_atoms = True
    if True:
        # try:

        # TODO --> need to loop over segment moltypes to check each one
        ###
        #moltype = allmoltype[i]
        moltype = m1.moltype()[0]  # for now ... early testing only       

        # if moltype.strip() not in ['protein','rna']:
        if(segbasis.lower() == 'backbone' or segbasis.lower() == 'heavy'):
            if(moltype == 'protein'):
                test_atoms = ["N", "CA", "C"]
            elif(moltype == 'rna'):
                test_atoms = ["P", "O5'", "C5'", "C3'", "O3'"]
            cutoff = 1.0

        elif(segbasis.lower() == 'all'):
            m1.set_average_vdw()
            atom_vdw = m1.atom_vdw()
            if any(None in sub_list for sub_list in atom_vdw):
                atom_list = ''
                name = m1.name()
                for i in xrange(m1.natoms()):
                    if(atom_vdw[i][0] == None):
                        atom_list += ' ' + name[i]
                error.append(
                    "the following atoms do not have vdW parameters; use different overlap basis : " + atom_list)  #NOT TESTED. There are no elements in the list that have "None" as a vdw parameter.
                return error
            else:
                check_atoms = False
                cutoff = 0.8
        else:
            string_basis = string.split(segbasis, ",")
            test_atoms = []
            print 'string_basis = ', string_basis
            # for i in xrange(len(string_basis)):
            for i in xrange(len(string_basis)):
                # for now ... early testing only
                test_atoms.append(string_basis[i].strip())
                cutoff = 2.0

# OPEN NEED OVERLAP ATOM CHECK FOR KEYWORDS HEAVY, BACKBONE, ETC.

        if check_atoms:
            for atom in test_atoms:
                if(atom not in name):
                    #error.append("overlap basis atom "+segbasis+" is not in your PDB file")
                    #error.append("overlap basis atom "+atom+" is not in your PDB file")
                    # return error
                    pass

    # TODO ... alter this to use the correct basis for each segment !!!

        check = 0
        coor = m1.coor()

        print 'checking overlap in initial structure using cutoff = ', cutoff
        print 'check = ', check
        check = overlap.overlap(coor[0], cutoff)
        print 'coor[0][0] = ', coor[0][0]
        print 'coor[0][-1] = ', coor[0][-1]
        print 'check = ', check

        if(check == 1):
            error.append(
                "Your original structure has overlapping atoms using your chosen basis.  If you fail to get accepted structures then either choose another basis or energy minimize your structure")
            print 'warning: ', error
            error = []

        # print m1.atom_vdw()

        # except:
        #	error.append("cannot parse overlap basis")
        #	return error


# TODO --> switching to non-single atom name overlap basis causes us to use hard-wired alignment basis in nmer_dihedral
###
#	for item in segbasis.split(','):
#		if item.strip() not in m1.name():
#			error.append('The segment alignment basis "'+item+'" is not found in the pdb file!')
#			return error

    for item in flpsegname.split(','):
        if item.strip() not in m1.segname():
            error.append('The flexible segment name "' + item +
                         '" is not found in the pdb file!')
            return error

    asegs = []
    segname = m1.segname()
    for s in range(len(segname)):
        tseg = segname[s]
        if(tseg not in asegs):
            asegs.append(tseg)
    numsegs = len(asegs)

    if(numsegs != nsegments):
        error.append(
            'the total number of segments entered does NOT match the number of segments in the pdb file')
        return error

    locvariables = ['resid']
    value, result = input_filter.get_pdb_stats(pdbfile, locvariables)
    resid = map(int, result[0])

#    if(resid[0] != 1):
#        error.append('amino acid residues in starting pdbfile '+pdbfile+' must start at resid = 1 : '+str(resid[0]))
#        return error


### HERE FOR SEGNAME RESID CHECK ####
### HERE FOR SEGNAME RESID CHECK ####
### HERE FOR SEGNAME RESID CHECK ####

    number_aa = resid[-1] - resid[0] + 1


###
# the following needs to loop over segment names
###
#    for i in xrange(resid[0],number_aa):
#        #ti=i+1
#        ti=i
#        if ti not in resid:
#            error.append('amino acid '+str(ti)+' is missing from pdbfile'+pdbfile)
#            return error

    lerror = check_segvariables(
        pdbfile, psegvariables, flpsegname, seglow, seghigh)
    if len(lerror):
        #error.append('failed to check segvariables in the complex filter!')
        error.append(lerror)
        return error

    return error
