import numpy as np
import os

def writeI_S(handles):
    # Check if the output file name ends with '.txt'
    if not handles['output_file_name'].endswith('.txt'):
        handles['output_file_name'] += '.txt'
    
    # Open file for writing
    with open(os.path.join(handles['output_file_path'], handles['output_file_name']), 'w') as file:
        # Write header information
        file.write('X-ray scattering calculations using labelled molecules\n')
        file.write('Intensity data\n')
        
        # Write parameters
        file.write(f"Number of discrete distances nD\t{handles['nDCalc']:.5g}\n")
        file.write(f"Maximum discrete distance DMax\t{handles['DMaxCalc']:.5g}\tAng\n")
        file.write(f"Maximum distance in molecule\t{handles['actDMax']:.5g}\tAng\n")
        file.write(f"Number of discrete S-values nS\t{handles['nSCalc']:.5g}\n")
        file.write(f"Maximum of S-value SMax\t{handles['SMaxCalc']:.5g}\t/Ang\n")
        file.write("Definition of S\tS=2*sin(theta)/lambda\n")
        
        # File names
        file.write(f"Input pdb file\t{handles['pdbFile']}\n")
        file.write(f"Scattering factor directory\t{handles['fDataPath']}\n")
        if handles['CRY']:
            file.write(f"CRYSOL intensity file\t{handles['CRYSOLFile']}\n")
        
        # Energies
        file.write('Energies / eV\n')
        for iEn in range(handles['nEnCalc']):
            file.write(f" \t{handles['EnCalc'][iEn]:.5g}\teV\n")
        file.write('\n')
        
        titles = [[], [], [], []]

        if handles['CRY']:
            titles[0].append('X-ray energy / eV')
            titles[1].append('Label i')
            titles[2].append('Label j')
            titles[3].append('S (/Ang)')
            for iEn in range(handles['nEnCalc']):
                titles[0].append(str(handles['EnCalc'][iEn]))
                titles[1].append('')
                titles[2].append('')
                titles[3].append(f"I(S) atoms at {str(handles['EnCalc'][iEn])} eV")

                titles[0].append('')
                titles[1].append('')
                titles[2].append('')
                titles[3].append('I(S) CRYSOL (in solvent)')
                
                titles[0].append('')
                titles[1].append('')
                titles[2].append('')
                titles[3].append('I(S) CRYSOL (in vacuum)')
            
                for ilblType in range(handles['nlblTypes']):
                    titles[0].append(str(handles['EnCalc'][iEn]))
                    titles[1].append(str(ilblType))
                    titles[2].append('')
                    titles[3].append(f"I(S) label {str(ilblType)} atoms at {str(handles['EnCalc'][iEn])} eV")
                
                for ilblType in range(handles['nlblTypes']):
                    for jlblType in range(ilblType, handles['nlblTypes']):
                        titles[0].append(str(handles['EnCalc'][iEn]))
                        titles[1].append(str(ilblType))
                        titles[2].append(str(jlblType))
                        titles[3].append(f"I(S) label {str(ilblType)}-{str(jlblType)} at {str(handles['EnCalc'][iEn])} eV")
                
                titles[0].append(str(handles['EnCalc'][iEn]))
                titles[1].append('')
                titles[2].append('')    
                titles[3].append(f"I(S) total (incl. atoms/labels) at {str(handles['EnCalc'][iEn])} eV")
                
                titles[0].append(str(handles['EnCalc'][iEn]))
                titles[1].append('')   
                titles[2].append('')
                titles[3].append(f"I(S) CRYSOL total in solvent at {str(handles['EnCalc'][iEn])} eV")

            formatStringTop = '\t'.join(['{}'] * len(titles[0])) + '\n'
            formatStringBot = '\t|'.join(['{}'] * len(titles[0])) + '\n'
            
            # Output the lines for this S-value
            file.write(formatStringTop.format(*titles[0]))
            file.write(formatStringTop.format(*titles[1]))
            file.write(formatStringTop.format(*titles[2]))
            file.write(formatStringBot.format(*titles[3]))

            for iS in range(handles['nSCalc'] + 1):
                line = []
                line.append(str(round(handles['SCalc'][iS], 4)))
                for iEn in range(handles['nEnCalc']):
                    line.append(str(format(float(handles['I_Satoms'][iEn, iS]), ".6f")))
                    toAdd = ''
                    if np.isnan(handles['I_SCRYatSolvS'][0][iS]):
                        toAdd = 'NaN'
                    else:
                        toAdd = str(format(float(handles['I_SCRYatSolvS'][0][iS]), ".6f"))
                    line.append(toAdd)

                    if np.isnan(handles['I_SCRYatVacS'][0][iS]):
                        toAdd = 'NaN'
                    else:
                        toAdd = str(format(float(handles['I_SCRYatVacS'][0][iS]), ".6f"))
                    line.append(toAdd)

                    for ilblType in range(handles['nlblTypes']):
                        line.append(str(format(float(handles['I_SlblAtoms'][iEn, ilblType, iS]), ".6f")))
                    for ilblType in range(handles['nlblTypes']):
                        for jlblType in range(ilblType, handles['nlblTypes']):
                            line.append(str(format(float(handles['I_SlblLbl'][iEn, ilblType, jlblType, iS]), ".6f")))
                    line.append(str(format(float(handles['I_Stot'][iEn, iS]), ".6f")))
                    if np.isnan(handles['I_SCRYtot'][iEn, iS]):
                        toAdd = 'NaN'
                    else:
                        toAdd = str(format(float(handles['I_SCRYtot'][iEn, iS]), ".6f"))
                    line.append(toAdd)

                formatString = '\t'.join(['{}'] * len(line)) + '\n'
                file.write(formatString.format(*line))
        else:
            titles[0].append('X-ray energy / eV')
            titles[1].append('Label i')
            titles[2].append('Label j')
            titles[3].append('S (/Ang)')
            for iEn in range(handles['nEnCalc']):
                titles[0].append(str(handles['EnCalc'][iEn]))
                titles[1].append('')
                titles[2].append('')
                titles[3].append(f"I(S) atoms at {str(handles['EnCalc'][iEn])} eV")
                for ilblType in range(handles['nlblTypes']):
                    titles[0].append(str(handles['EnCalc'][iEn]))
                    titles[1].append(str(ilblType))
                    titles[2].append('')
                    titles[3].append(f"I(S) label {str(ilblType)} /atoms at {str(handles['EnCalc'][iEn])} eV")
                for ilblType in range(handles['nlblTypes']):
                    for jlblType in range(ilblType, handles['nlblTypes']):
                        titles[0].append(str(handles['EnCalc'][iEn]))
                        titles[1].append(str(ilblType))
                        titles[2].append(str(jlblType))
                        titles[3].append(f"I(S) label {str(ilblType)}-{str(jlblType)} at {str(handles['EnCalc'][iEn])} eV")
                titles[0].append(str(handles['EnCalc'][iEn]))
                titles[1].append('')
                titles[2].append('')
                titles[3].append(f"I(S) total (incl. atoms/labels) at {str(handles['EnCalc'][iEn])} eV")
            
            formatStringTop = '\t'.join(['{}'] * len(titles[0])) + '\n'
            formatStringBot = '\t|'.join(['{}'] * len(titles[0])) + '\n'
            
            # Output the lines for this S-value
            file.write(formatStringTop.format(*titles[0]))
            file.write(formatStringTop.format(*titles[1]))
            file.write(formatStringTop.format(*titles[2]))
            file.write(formatStringBot.format(*titles[3]))

            for iS in range(handles['nSCalc'] + 1):
                line = []
                line.append(str(round(handles['SCalc'][iS], 4)))
                for iEn in range(handles['nEnCalc']):
                    line.append(str(format(int(handles['I_Satoms'][iEn, iS]), ".6f")))
                    for ilblType in range(handles['nlblTypes']):
                        line.append(str(format(int(handles['I_SlblAtoms'][iEn, ilblType, iS]), ".6f")))
                    for ilblType in range(handles['nlblTypes']):
                        for jlblType in range(ilblType, handles['nlblTypes']):
                            line.append(str(format(int(handles['I_SlblLbl'][iEn, ilblType, jlblType, iS]), ".6f")))
                    line.append(str(format(int(handles['I_Stot'][iEn, iS]), ".6f")))

                formatString = '\t'.join(['{}'] * len(line)) + '\n'
                file.write(formatString.format(*line))
    return

def write_used_vars(handles):
    try:
        if not handles['output_vars_file_name'].endswith('.txt'):
            handles['output_vars_file_name'] += '.txt'

        file_path = os.path.join(handles['output_vars_file_path'], handles['output_vars_file_name'])

        with open(file_path, 'w') as file:
            # Prints out truth statements and energy
            file.write(f"CRY={str(handles['CRY'])}\n")
            file.write(f"lblXal={str(handles['lblXal'])}\n")
            file.write(f"fPairsCalc={str(handles['fPairsCalc'])}\n")
            file.write(f"I_StotDone={str(handles['I_StotDone'])}\n")
            file.write(f"I_SCRYtotDone={str(handles['I_SCRYtotDone'])}\n")

            # Prints out parameters
            file.write(f"DMax={handles['DMax']}\n")
            file.write(f"nD={handles['nD']}\n")
            file.write(f"SMax={handles['SMax']:.4f}\n")
            file.write(f"nS={handles['nS']}\n")
            file.write(f"err={handles['err']:.4f}\n")
            file.write(f"DMaxInv={handles['DMaxInv']}\n")
            file.write(f"nDInv={handles['nDInv']}\n")
            file.write(f"actDMax={handles['actDMax']:.4f}\n")

            # Prints out label and atom information
            output_line = ' '.join(map(str, handles['lblType']))
            file.write(f"lblType= {output_line}\n")

            output_line = ' '.join(map(str, handles['numlblType']))
            file.write(f"numlblType= {output_line}\n")

            output_line = ' '.join(map(str, handles['lblRad']))
            file.write(f"lblRad= {output_line}\n")

            output_line = ' '.join(map(str, handles['lblNum']))
            file.write(f"lblNum= {output_line}\n")

            file.write('tbl_numType=\n')
            for atom in range(len(handles['atomSet'])):
                file.write(f"{handles['tbl_numType'][atom][0]} {handles['tbl_numType'][atom][1]}\n")

            # Prints out EnCalc
            output_line = ' '.join(map(str, handles['EnCalc']))
            file.write(f"EnCalc=\n{output_line}\n")

            # Outputs results of pairScattFacs
            if handles['fPairsCalc']:
                # Outputs f
                file.write('fCalc=\n')
                for i_type in range(len(handles['fCalc'][0])):
                    file.write(f"{handles['atomSet'][i_type]}\n")
                    for i_en in range(len(handles['fCalc'])):
                        f_real = handles['fCalc'][i_en][i_type].real
                        f_imag = handles['fCalc'][i_en][i_type].imag
                        file.write(f"{f_real}+{f_imag}j\n")

                # Outputs fPair
                file.write('fPairCalc=\n')
                for i_type in range(len(handles['fPairCalc'][0])):
                    for j_type in range(len(handles['fPairCalc'][0][0])):
                        file.write(f"{handles['atomSet'][i_type]}-{handles['atomSet'][j_type]}\n")
                        output_line = ' '.join(map(str, handles['fPairCalc'][:, i_type, j_type]))
                        file.write(f"{output_line}\n")

                # Outputs fMeanAtoms
                file.write('meanfAtCalc=\n')
                for i in range(len(handles['meanfAtCalc'][0])):
                    output_line = ' '.join(map(str, handles['meanfAtCalc'][:, i]))
                    file.write(f"{output_line}\n")

            # Prints out SCalc
            file.write('SCalc\n')
            for i_s in range(handles['nSCalc'] + 1):
                file.write(f"{handles['SCalc'][i_s]:12.5g}\n")

            # Prints out I_Stot
            file.write('I_Stot\n')
            for i_s in range(handles['nSCalc'] + 1):
                output_line = ' '.join(map(str, handles['I_Stot'][:, i_s]))
                file.write(f"{output_line}\n")

            # Prints out I_SCRYtot
            if handles['CRY']:
                file.write('I_SCRYtot\n')
                for i_s in range(handles['nSCalc'] + 1):
                    output_line = ' '.join(map(str, handles['I_SCRYtot'][:, i_s]))
                    file.write(f"{output_line}\n")

    except IOError:
        print('Failed to open file - is it open already?')
    return