
import numpy
import sasmol as system
import math
import sys
import os
import scipy
from scipy.constants import pi
from scipy.special import sinc
from scipy.interpolate import interp1d


def readfData(handles):
    # root = tk.Tk()
    # root.withdraw()  # Hide the main window
    # fDataPath = filedialog.askdirectory(initialdir=os.getcwd(), title='Identify the location of the scattering amplitude text files')


    fidfData = []
    for atomType in handles['atomSet']:
        try:
            f = open(os.path.join(handles['fDataPath'], f'{atomType}.txt'), 'r')
            fidfData.append(f)
        except IOError:
            print('Cannot open f-data files *.txt')
            return False, handles
    # No file errors, read in data
    handles['fData'] = [[None, None, None] for _ in range(len(fidfData))]
    errorlg = False

    for iFile, f in enumerate(fidfData):
        next(f)  # Discard the first line of titles
        data = numpy.loadtxt(f, dtype=float, unpack=True)
        if numpy.any(numpy.isnan(data)):
            errorlg = True
            break
        handles['fData'][iFile] = data

    if errorlg:
        print("Error: NaN values found in data")

    for f in fidfData:  # Close all files
        f.close()

    return handles

def pairScattFacs(energy, handles):
    nEn = len(energy)
    nTypes = len(handles['fData'])
    f = numpy.zeros((nEn, nTypes), dtype=complex)
    fPair = numpy.zeros((nEn, nTypes, nTypes))
    meanfAtoms = numpy.zeros((nEn, 2))

    for iEn in range(nEn):
        for iType in range(nTypes):
            # real_interp = interp1d(handles['fData'][iType][0], handles['fData'][iType][1])
            # imag_interp = interp1d(handles['fData'][iType][0], handles['fData'][iType][2])
            # f[iEn, iType] = complex(real_interp(energy[iEn]), imag_interp(energy[iEn]))
            real_interp_value = numpy.interp(energy[iEn], handles['fData'][iType][0], handles['fData'][iType][1])
            imag_interp_value = numpy.interp(energy[iEn], handles['fData'][iType][0], handles['fData'][iType][2])
            f[iEn, iType] = complex(real_interp_value, imag_interp_value)

        for iType in range(nTypes):
            for jType in range(nTypes):
                fPair[iEn, iType, jType] = round(2 * numpy.real(f[iEn, iType] * numpy.conj(f[iEn, jType])), 6)

    numType = handles['tbl_numType']
    numEstSet = True
    numSet = True
    if numType[0][2] is None:
        numEstSet = False
        meanfAtoms[:, 1] = 1.0
    if numType[0][1] is None:
        numSet = False
        meanfAtoms[:, 0] = 1.0

    for iEn in range(nEn):
        for iType in range(len(handles['atomSet'])):
            if numSet:
                meanfAtoms[iEn, 0] += numType[iType][1] * abs(f[iEn, iType])**2
            if numEstSet:
                meanfAtoms[iEn, 1] += numType[iType][2] * abs(f[iEn, iType])**2

        for iEn in range(1, nEn):
            meanfAtoms[iEn, 0] /= meanfAtoms[0, 0]
            meanfAtoms[iEn, 1] /= meanfAtoms[0, 1]
    meanfAtoms[0, :] = 1.0

    handles['fPairsCalc'] = True

    #output
    # output_file_path = '/Users/byc1/Documents/GitHub/asaxs_code/src/scat_label_and_sum/python_prototype/'
    # output_file_name = f"fMeanAtomsTest"


    # if not output_file_name.endswith('.txt'):
    #     output_file_name += '.txt'

    # # For f
    # with open(os.path.join(output_file_path, output_file_name), 'w') as file:
    #     for iEn in range(len(f)):
    #         for iType in range(len(f[iEn])):
    #             file.write(f"f[{iEn}, {iType}] = {f[iEn, iType]}\n")


    # # For fPair
    # with open(os.path.join(output_file_path, output_file_name), 'w') as file:
    #     for iEn in range(fPair.shape[0]):  # Energy index
    #         for iType in range(fPair.shape[1]):  # Type index for the first atom
    #             for jType in range(fPair.shape[2]):  # Type index for the second atom
    #                 # Step 4: Write the value to the file
    #                 file.write(f"fPair[{iEn}, {iType}, {jType}] = {fPair[iEn, iType, jType]}\n")

    # # For meanfAtoms
    # with open(os.path.join(output_file_path, output_file_name), 'w') as file:
    #     for iEn in range(meanfAtoms.shape[0]):
    #         file.write(f"meanfAtoms[{iEn}, {0}] = {meanfAtoms[iEn, 0]}, meanfAtoms[{iEn}, {1} = {meanfAtoms[iEn, 1]}\n")

    return f, fPair, meanfAtoms

## GITHUB COPILOT CONVERSION (JEC):

#import numpy as np
#from scipy.interpolate import interp1d

def pair_scatt_facs(energy, handles):
    """
    Create an array of 2*Re(fi*conj(fj) and fi for all possible atom types, atom types pairs, and energies.
    Input 'energy' is an array of energy values in keV.
    Uses the fData array of scattering factors held in the handles structure.
    Outputs fPair(energy, iType, jType) for each pair of atoms and energy,
    and f(iEn, iType) and meanfAtoms(iEn, 2) (an average atom pair-scattering factor, 1:pdb data, 2:est).
    """
    nEn = len(energy)
    nTypes = len(handles['fData'])
    f = numpy.zeros((nEn, nTypes), dtype=complex)
    fPair = numpy.zeros((nEn, nTypes, nTypes))

    for iEn in range(nEn):
        for iType in range(nTypes):  # Calculate the f at the required energy for each type
            real_interp = interp1d(handles['fData'][iType][0], handles['fData'][iType][1], kind='linear', bounds_error=True)
            imag_interp = interp1d(handles['fData'][iType][0], handles['fData'][iType][2], kind='linear', bounds_error=True)
            f[iEn, iType] = complex(real_interp(energy[iEn]), imag_interp(energy[iEn]))

        # Now calculate the combined scattering factors of all possible pairs
        for iType in range(nTypes):
            for jType in range(nTypes):
                fPair[iEn, iType, jType] = 2 * numpy.real(f[iEn, iType] * numpy.conj(f[iEn, jType]))

    # Calculate the mean abs(f)*abs(f) for self-scattering by all types of atom in molecule weighted by numbers of atoms of that type
    meanfAtoms = numpy.zeros((nEn, 2))
    numType = handles['tbl_numType']
    numEstSet = True
    numSet = True

    if not numType[0][2]:
        numEstSet = False
        meanfAtoms[:, 1] = 1.0

    if not numType[0][1]:
        numSet = False
        meanfAtoms[:, 0] = 1.0

    for iEn in range(nEn):
        for iType in range(len(handles['atomSet'])):
            if numSet:
                meanfAtoms[iEn, 0] += numType[iType][1] * abs(f[iEn, iType]) * abs(f[iEn, iType])
                #pprint.pp(numType[iType][1])

            elif numEstSet:
                meanfAtoms[iEn, 1] += numType[iType][2] * abs(f[iEn, iType]) * abs(f[iEn, iType])
                #pprint(numType[iType][1])

            # meanfAtoms[iEn, 1] is averaged over the counted number of each type
            # meanfAtoms[iEn, 2] is averaged over the estimated number of each type (user input)

    # Ratio each value by the value at the first energy
    for iEn in range(1, nEn):
        meanfAtoms[iEn, 0] = meanfAtoms[iEn, 0] / meanfAtoms[0, 0]
        meanfAtoms[iEn, 1] = meanfAtoms[iEn, 1] / meanfAtoms[0, 1]

    meanfAtoms[0, :] = 1.0
    fPairsCalc = True

    return f, fPair, meanfAtoms, fPairsCalc  # Adjust the return values as needed


def calcI_Satoms(handles):
    errorlg = False
    nEn = handles['nEnCalc']
    nD = handles['nDCalc']
    nS = handles['nSCalc']

    deltaD = handles['deltaD']
    SCalc = handles['SCalc']
    actDMax = handles['actDMax']
    fPairs = handles['fPairCalc']
    nAtoms = len(handles['atomData'])

    amplDatoms = numpy.zeros((nEn, nD + 1))  # all atoms with each other (incl D=0 -> nD+1)
    handles['I_Satoms'] = numpy.zeros((nEn, nS + 1))  # all atoms with each other

    for iAtom in range(nAtoms):  # for each atom
        iType = handles['atomData'][iAtom][1]
        coordsi = handles['atomData'][iAtom][3]
        xi, yi, zi = coordsi

        for iEn in range(nEn):  # assign to 1st D-bin for each energy
            amplDatoms[iEn, 0] += fPairs[iEn, iType, iType]

        # For testing purposes
        # print(f"amplDatoms: {amplDatoms[:,0]}")

        for jAtom in range(iAtom + 1, nAtoms):  # count pairs EXcluding self-scatter
            coordsj = handles['atomData'][jAtom][3]
            jType = handles['atomData'][jAtom][1]

            D = math.sqrt((coordsj[0] - xi)**2 + (coordsj[1] - yi)**2 + (coordsj[2] - zi)**2)
            actDMax = max(actDMax, D)

            iD = math.ceil(D / deltaD)  # the element number for the D bin

            if iD > nD:
                errorlg = True
                raise Exception('DMax needs to be increased')
                # Replace with appropriate error handling
            else:
                for iEn in range(nEn):  # assign to D-bin for each energy
                    amplDatoms[iEn, iD+1] += fPairs[iEn, iType, jType]

            if errorlg:
                break

        if errorlg:
            break

    if not errorlg:
        for iEn in range(nEn):
            for iD in range(nD + 1):
                SD = 2 * handles['SCalc'] * iD * handles['deltaD']
                handles['I_Satoms'][iEn,:] += amplDatoms[iEn, iD] * numpy.sinc(SD)

    handles['actDMax'] = actDMax
    return handles, errorlg

def f_Ssphere(S, R, n):
    """
    Calculate the S-dependent factor to multiply the atom scattering factor
    for a finite-size sphere.
    Inputs:
    - S: array of S-values
    - R: radius of nanocrystal
    - n: number of atoms in nanocrystal
    """
    nBins = 1
    RSet = numpy.array([R])
    f_fac = numpy.zeros_like(S)

    for iBin in range(nBins):
        x = 2 * numpy.pi * RSet[iBin] * S
        f_fac = n * numpy.ones_like(S)
        lnonZero = numpy.where(x > 1e-10)
        f_fac[lnonZero] = 3 * n * (numpy.sin(x[lnonZero]) - x[lnonZero] * numpy.cos(x[lnonZero])) / (x[lnonZero] ** 3)

    return f_fac


def calcI_SlblAtoms(handles):
    errorlg = False
    nEn = handles['nEnCalc']
    nD = handles['nDCalc']
    nS = handles['nSCalc']
    nlblTypes = handles['nlblTypes']
    nLabels = handles['nLabels']

    amplDlblAt = numpy.zeros((nEn, nlblTypes, nD+1))
    handles['I_SlblAtoms'] = numpy.zeros((nEn, nlblTypes, nS+1))

    actDMax = handles['actDMax']
    fPairs = handles['fPairCalc']
    deltaD = handles['deltaD']
    SCalc = handles['SCalc']

    if not errorlg:
        for iLabel in range(nLabels):
            ilblCount = handles['labelData'][iLabel][4]
            iType = handles['labelData'][iLabel][1]
            coordsi = handles['labelData'][iLabel][3]
            xi, yi, zi = coordsi
            for iEn in range(nEn):
                amplDlblAt[iEn, ilblCount, 0] += fPairs[iEn, iType, iType]

            for jAtom in range(handles['nAtoms']):
                coordsj = handles['atomData'][jAtom][3]
                jType = handles['atomData'][jAtom][1]
                D = numpy.sqrt((coordsj[0] - xi)**2 + (coordsj[1] - yi)**2 + (coordsj[2] - zi)**2)
                actDMax = max(actDMax, D)
                iD = math.ceil(D / deltaD)
                if iD > nD:
                    errorlg = True
                    print('DMax needs to be increased')
                    break
                else:
                    for iEn in range(nEn):
                        amplDlblAt[iEn, ilblCount, iD] += fPairs[iEn, iType, jType]

    handles['actDMax'] = actDMax

    f_shape = numpy.ones((nlblTypes, len(SCalc)))
    #print(numpy.shape(f_shape))
    if handles['lblXal']:
        for ilbl in range(nlblTypes):
            f_shape[ilbl, :] = f_Ssphere(SCalc, handles['lblRadCalc'][ilbl], handles['lblNumCalc'][ilbl])

    if not errorlg:
        for iEn in range(nEn):
            for ilblCount in range(nlblTypes):
                for iD in range(nD+1):
                    SD = 2 * SCalc * iD * deltaD
                    # did not know how to implement the reshape method, works without it but will result in errors with it. Shape may already match?
                    handles['I_SlblAtoms'][iEn, ilblCount, :] += amplDlblAt[iEn, ilblCount, iD] * sinc(SD) * f_shape[ilblCount, :]

    return handles, errorlg

def distPairs(nEn, iType, jType, coordsi, coordsj, fPair, amplD, handles):
    """
    Assign the amplitude scattering array as function of distance.
    Calculates the separation of two atoms/labels and assigns their scattering
    amplitude to the appropriate distance in the array.

    :param nEn: The energy element.
    :param iType: The type of atom or label for the i-th atom/label.
    :param jType: The type of atom or label for the j-th atom/label.
    :param coordsi: The coordinates of the i-th atom/label.
    :param coordsj: The coordinates of the j-th atom/label.
    :param fPair: The array containing the pair scattering factors.
    :param amplD: The array to which the scattering should be assigned.
    :param handles: The handles array holding actDMax (tracking maximum separation) and D values.
    :return: Updated handles, amplD, and errorlg indicating if an error occurred.
    """
    errorlg = False
    # Calculate distance between coordsi and coordsj
    D = numpy.sqrt((coordsj[0] - coordsi[0])**2 + (coordsj[1] - coordsi[1])**2 + (coordsj[2] - coordsi[2])**2)
    handles['actDMax'] = max(handles['actDMax'], D)  # Save the biggest value

    if handles['actDMax'] > handles['DMaxCalc']:
        errorlg = True
        # In actual implementation, display error dialog here
        print('DMax needs to be increased')  # Placeholder for error dialog
    else:
        if D < 1e-5:
            iD = 0  # self-scattering
        else:
            iD = int(numpy.ceil(D / handles['deltaD']))  # the element number for the D bin

        for iEn in range(nEn):  # Assign to D-bin for each energy
            amplD[iEn, iD] += fPair[iEn, iType, jType]  # Assign the scattering to the relevant D bin

    return handles, amplD, errorlg



def calcI_SlblLbl(handles):
    errorlg = False
    nEn = handles['nEnCalc']
    nD = handles['nDCalc']
    nS = handles['nSCalc']
    nlblTypes = handles['nlblTypes']
    nLabels = handles['nLabels']

    amplDlbl = numpy.zeros((nEn, nlblTypes, nlblTypes, nD + 1))  # each label with another label
    handles['I_SlblLbl'] = numpy.zeros((nEn, nlblTypes, nlblTypes, nS + 1))

    if not errorlg:
        for iLabel in range(nLabels):
            ilblCount = handles['labelData'][iLabel][4]  # Adjusted for 0-based indexing
            for jLabel in range(iLabel, nLabels):
                jlblCount = handles['labelData'][jLabel][4]  # Adjusted for 0-based indexing
                ordLblCount = sorted([ilblCount, jlblCount])
                ilbl, jlbl = ordLblCount
                handles, amplDlbl[:, ilbl, jlbl, :], errorlg = distPairs(nEn, handles['labelData'][iLabel][1], handles['labelData'][jLabel][1], handles['labelData'][iLabel][3], handles['labelData'][jLabel][3], handles['fPairCalc'], amplDlbl[:, ilbl, jlbl, :], handles)
                if errorlg:
                    break

    f_shape = numpy.ones((nlblTypes, len(handles['SCalc'])))
    if handles['lblXal']:
        for ilbl in range(nlblTypes):
            f_shape[ilbl, :] = f_Ssphere(handles['SCalc'], handles['lblRadCalc'][ilbl], handles['lblNumCalc'][ilbl])

    if not errorlg:
        for iEn in range(nEn):
            for ilblCount in range(nlblTypes):
                for jlblCount in range(ilblCount, nlblTypes):
                    for iD in range(nD + 1):
                        SD = 2 * handles['SCalc'] * iD * handles['deltaD']
                        handles['I_SlblLbl'][iEn, ilblCount, jlblCount, :] += amplDlbl[iEn, ilblCount, jlblCount, iD] * numpy.sinc(SD) * f_shape[ilblCount, :]

    return handles, errorlg

def calcI_Stotal(handles):
    try:
        errorlg = False
        nEn = handles['nEnCalc']
        nS = handles['nSCalc']
        handles['I_Stot'] = numpy.zeros((nEn, nS + 1))
        handles['I_SCRYtot'] = numpy.zeros((nEn, nS + 1))
        handles['I_SCRYtotVac'] = numpy.zeros((nEn, nS + 1))
        sizeS = handles['I_Stot'].shape[0]
        nlblTypes = handles['nlblTypes']
        for iEn in range(nEn):
            handles['I_Stot'][iEn, :] = handles['I_Satoms'][iEn, :]  # start with the atoms only
            if 'CRY' in handles and handles['CRY']:  # atoms +solvent from CRYSOL
                handles['I_SCRYtot'][iEn, :] = handles['I_SCRYatSolvS']
                handles['I_SCRYtotVac'][iEn, :] = handles['I_SCRYatVacS']
            for ilblType in range(nlblTypes):
                if 'CRY' in handles and handles['CRY']:  # only calc if CRYSOL data being used
                    handles['I_SCRYtot'][iEn, :] += handles['I_SlblAtoms'][iEn, ilblType, :]  # label-atom scattering with CRYSOL
                    handles['I_SCRYtotVac'][iEn, :] += handles['I_SlblAtoms'][iEn, ilblType, :]  # label-atom scattering with CRYSOL
                handles['I_Stot'][iEn, :] += handles['I_SlblAtoms'][iEn, ilblType, :]  # label-atom scattering with atom scattering
                for jlblType in range(ilblType, nlblTypes):  # label-label scattering
                    if 'CRY' in handles and handles['CRY']:
                        handles['I_SCRYtot'][iEn, :] += handles['I_SlblLbl'][iEn, ilblType, jlblType, :]  # label-label scattering with CRYSOL
                        handles['I_SCRYtotVac'][iEn, :] += handles['I_SlblLbl'][iEn, ilblType, jlblType, :]  # label-label scattering with CRYSOL
                    handles['I_Stot'][iEn, :] += handles['I_SlblLbl'][iEn, ilblType, jlblType, :]  # label-label scattering with atom scattering
    except Exception as MError1:
        errorlg = True

    return handles, errorlg


def calculate_scattering_and_sum(handles):

###  original python prototype by Benson Chan

    handles = input_data.readPDB(handles)

    if crysol_file_flag:
        handles = input_data.read.CRYSOL(handles)

    handles = readfData(handles)

    errorlg = False

    # Save the parameters as set for this calculation
    # because they may later be changed by the user in the GUI

    handles['DMaxCalc'] = handles['DMax']
    handles['nDCalc'] = handles['nD']
    handles['EnCalc'] = handles['XrayEn']
    handles['nEnCalc'] = len(handles['XrayEn'])
    handles['SMaxCalc'] = handles['SMax']
    handles['nSCalc'] = handles['nS']
    nEn = handles['nEnCalc']  # Number of different energies to be used
   

    # There are handles['nD']+1 distance bins, including D=0 and D=handles['DMax']
    # All distances in Angstroms
    handles['deltaD'] = handles['DMaxCalc'] / handles['nDCalc']  # in Angstroms
    deltaS = handles['SMaxCalc'] / handles['nSCalc']
    handles['SCalc'] = numpy.linspace(0, handles['SMaxCalc'], handles['nSCalc'] + 1)

    # Save the nanocrystal settings
    handles['lblRadCalc'] = handles['lblRad']  # Save the radius and number atoms settings
    handles['lblNumCalc'] = handles['lblNum']

    # Calculate combined scattering factor
#    handles['fCalc'], handles['fPairCalc'], handles['meanfAtCalc'] = pairScattFacs(handles['EnCalc'], handles)

    handles['fCalc'], handles['fPairCalc'], handles['meanfAtCalc'], handles['fPairsCalc'] = pair_scatt_facs(handles['EnCalc'], handles)

    import pprint

    print('fCalc')

    pprint.pp(handles['fCalc'])

    print('fPairCalc')

    pprint.pp(handles['fPairCalc'])

    print('meanfAtCalc')

    pprint.pp(handles['meanfAtCalc'])

    # Calculate the atom-atom scattering I(S)
    handles, errorlg = calcI_Satoms(handles)
    if errorlg:
        raise Exception('Atom scatter intensity calc error')

    #print('I_Satoms')

    #pprint.pp(handles['I_Satoms'][-1:])

    #print('DONE')

    # Calculate the label-atom scattering I(S)
    handles, errorlg = calcI_SlblAtoms(handles)
    if errorlg:
        raise Exception('Label-scatter intensity calc error')


    # Calculate the label-label scattering I(S)
    handles, errorlg = calcI_SlblLbl(handles)
    if errorlg:
        raise Exception('Label-label intensity calc error')

    handles, errorlg = calcI_Stotal(handles)
    if errorlg:
        raise Exception('Total scatter intensity calc error')

    handles['I_StotDone'] = True
    handles['I_SCRYtotDone'] = True

    """
    
    # WRITEI_S METHOD WORKS ONLY FOR THE CASE WHEN CRYSOL DATA IS USED, ADD ADDITIONAL CASE WHEN NECESSARY

    #if handles['chk_addErr']:
        #handles['I_Stot'], handles['I_SCRYtot']=addErrors(handles['I_Stot'], handles['I_SCRYtot'], handles['err'], handles['CRY'])
    #writeI_S(handles)
    #print('Finished')

    """

    return



if __name__ == "__main__":

    import sassie.calculate.asaxs.asaxs_methods.prototype_testing_data.scattering_label_and_sum_data.input_data as input_data

    crysol_file_flag = False

    handles = input_data.initialize_input_data()

    calculate_scattering_and_sum(handles)

    print('done')

    #print('handles = ', handles)

    pass

