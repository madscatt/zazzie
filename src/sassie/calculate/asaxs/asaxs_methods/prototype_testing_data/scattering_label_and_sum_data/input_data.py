import numpy
import os

# Developer path
developer_path = os.path.join('/Users', 'byc1', 'Documents', 'GitHub')
developer_path = os.path.join('/Users', 'curtisj', 'git_working_copies')

print("developer path = ", developer_path)

def initialize_input_data():
    # Parameters to match the example data
    atomSet=['C', 'N', 'O', 'S', 'P', 'H', 'ZN', 'AU', 'AN', 'TB','PT']
    DMax=60
    nD=60
    SMax=0.1
    nS=1000
    err=0.01
    DMaxInv=180
    nDInv=90
    XrayEn = [11600, 11800, 12000, 12200, 12400]
    listEnergies = numpy.zeros(len(XrayEn))
    for iEn in range(len(XrayEn)):
        listEnergies[iEn] = str(XrayEn[iEn])
    actDMax=0
    fPairsCalc=False
    I_StotDone=False
    I_SCRYtotDone=False
    CRY=False
    I_SExpttotDone=False

    P_DCalcDone = False
    P_DCRYDone = False
    P_DExptDone = False
    err = 0.01

    CRYSOLFile = '*.int'
    pdbFile = '*.pdb'
    #fDataPath = os.getcwd()
    fDataPath = '' #CHANGE LATER TO IN THE CASE OF ERROR OCCURING
    handles = {'atomSet': atomSet,
           'DMax': DMax,
           'nD': nD,
           'SMax': SMax,
           'nS': nS,
           'err': err,
           'DMaxInv': DMaxInv,
           'nDInv': nDInv,
           'XrayEn': XrayEn,
           'listEnergies': listEnergies,
           'actDMax': actDMax,
           'fPairsCalc': fPairsCalc,
           'I_StotDone': I_StotDone,
           'I_SCRYtotDone': I_SCRYtotDone,
           'CRY': CRY,
           'I_SExpttotDone': I_SExpttotDone,
           'P_DCalcDone': P_DCalcDone,
           'P_DCRYDone': P_DCRYDone,
           'P_DExptDone': P_DExptDone,
           'CRYSOLFile': CRYSOLFile,
           'pdbFile': pdbFile,
           'fDataPath': fDataPath,
           'err': err}


    handles['pdbFilePath'] = developer_path + '/asaxs_code/src/scat_label_and_sum/python_prototype/test_pdb_inputs/'
    handles['pdbFile'] = 'python_output_read_pdb_DNA10Au2.txt'
    handles['CRYSOLFile'] = developer_path + '/asaxs_code/src/scat_label_and_sum/python_prototype/test_CRYSOL_inputs/matlab_output_read_crysol_DNA10Au2.txt'
    handles['fDataPath'] = developer_path + '/asaxs_code/src/scat_label_and_sum/matlab/scattering_factors/'
    handles['output_file_path'] = developer_path + '/asaxs_code/src/scat_label_and_sum/python_prototype/python_prototype_test_scat_label_sum_outputs/'
    handles['output_file_name'] = f"python_prototype_output_calc_tot_I_{handles['pdbFile'].replace('python_output_read_pdb_', '')}"

    handles['setRad'] = 7
    handles['setNum'] = 78

    # HARDCODED, SHOULD BE OFFERED AS OPTIONS IN GUI
    handles['lblXal'] = True
    handles['CRY'] = False
    handles['chk_addErr'] = False

    return handles

def readPDB(handles):
    # root = tk.Tk()
    # root.withdraw()
    # handles['pdbFile'] = filedialog.askopenfilename(initialdir = os.getcwd(), title = 'Select PDB file', filetypes = (('txt files', '*.txt'), ('all files', '*.*')))
    # if handles['pdbFile'] == '':
    #     raise Exception('No file selected')
    # BENSON, HARDCODED FOR TESTING

    handles['atomData'] = []
    handles['labelData'] = []
    handles['nAtomsType'] = numpy.zeros(len(handles['atomSet']))
    handles['nLabels'] = 0
    handles['lblType'] = []
    handles['nAtoms'] = 0
    handles['numlblType'] = []

    with open((handles['pdbFilePath'] + handles['pdbFile']), 'r') as fid:
        next(fid)
        for line in fid:
            if line != 'Label Data:\n':
                parts = line.split()
                if len(parts) == 6:
                    handles['atomData'].append([parts[0], int(parts[1])-1, parts[2], [float(parts[3]), float(parts[4]), float(parts[5])]])
                    handles['nAtomsType'][int(parts[1])-1] += 1
                else:
                    raise Exception('Invalid format')
            else:
                break

        iLabel = 0
        for line in fid:
            parts = line.split()
            if len(parts) == 6:
                handles['labelData'].append([parts[0], int(parts[1])-1, parts[2], [float(parts[3]), float(parts[4]), float(parts[5])]])
                handles['nLabels'] += 1
                if handles['labelData'][iLabel][1] not in handles['lblType']:
                    handles['lblType'].append(handles['labelData'][iLabel][1])
                    handles['numlblType'].append(1)
                else:
                    handles['numlblType'][handles['lblType'].index(handles['labelData'][iLabel][1])] += 1
                handles['labelData'][iLabel].append(handles['lblType'].index(handles['labelData'][iLabel][1]))
                iLabel = iLabel + 1
            else:
                raise Exception('Invalid format')

    handles['nlblTypes'] = len(handles['lblType'])
    handles['lblRad'] = handles['setRad']*numpy.ones(numpy.shape(handles['lblType'])) # HARDCODED, SHOULD BE OFFERED AS OPTIONS IN GUI
    handles['lblNum'] = handles['setNum']*numpy.ones(numpy.shape(handles['lblType'])) # HARDCODED, SHOULD BE OFFERED AS OPTIONS IN GUI

    # print(handles['lblType'])
    # print(handles['lblRad'])
    # print(handles['lblNum'])

    nTypes = len(handles['atomSet'])
    numTypeTbl = []  # Create the data for the table

    for iType in range(nTypes):
        atomType = handles['atomSet'][iType]
        atomCount = int(handles['nAtomsType'][iType])  # The counted number of atoms of each type from pdb file
        numTypeTbl.append([atomType, atomCount, None])

    handles['tbl_numType'] = numTypeTbl

    return handles

def readCRYSOL(handles):
    # root = tk.Tk()
    # root.withdraw()
    # handles['CRYSOLFile'] = filedialog.askopenfilename(initialdir = os.getcwd(), title = 'Select CRYSOL file', filetypes = (('txt files', '*.txt'), ('all files', '*.*')))
    # BENSON, HARDCODED FOR TESTING

    data = []
    with open(handles['CRYSOLFile'], 'r') as fid:
        for _ in range(3):
            next(fid)
        for line in fid:
            parts = line.split()
            if len(parts) == 3:
                data.append([float(parts[0]), float(parts[1]), float(parts[2])])
            else:
                raise Exception('Invalid format')
    nS = len(data)
    handles['CRYSOL_S'] = [row[0] for row in data]
    handles['I_SCRYatSolv'] = np.array([row[1] for row in data]).reshape(1, nS)
    handles['I_SCRYatVac'] = np.array([row[2] for row in data]).reshape(1, nS)
    # handles['nS'] = nS







    return handles



