import numpy

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

    return handles
