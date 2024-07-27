def write_input_data(handles):
    with open('input_data_calcI_Satoms_python.txt', 'w') as file:
        # Write nEnCalc
        file.write(f'nEnCalc: {handles["nEnCalc"]}\n')
        
        # Write nDCalc
        file.write(f'nDCalc: {handles["nDCalc"]}\n')
        
        # Write nSCalc
        file.write(f'nSCalc: {handles["nSCalc"]}\n')
        
        # Write deltaD
        file.write(f'deltaD: {handles["deltaD"]}\n')
        
        # Write SCalc
        SCalc = handles['SCalc']
        for element in SCalc:
            file.write(f'SCalc: {element:.6f}\n')

        # Write actDMax
        file.write(f'actDMax: {handles["actDMax"]}\n')
        
        # Write fPairCalc
        file.write('fPairCalc:\n')
        fPairCalc = handles['fPairCalc']
        for i in range(fPairCalc.shape[0]):
            for j in range(fPairCalc.shape[1]):
                for k in range(fPairCalc.shape[2]):
                    file.write(f'{fPairCalc[i, j, k]:.6f} ')
                file.write('\n')
        
        # Write nAtoms
        nAtoms = len(handles['atomData'])
        file.write(f'nAtoms: {nAtoms}\n')
        file.write(f'atomData:\n')
        
        for i in range(len(handles['atomData'])):
            file.write(f'Atom {i + 1}:\n')  # MATLAB indexing starts at 1
            file.write(f'Type: {handles["atomData"][i][1] + 1}\n')
            coordinates = handles['atomData'][i][3]
            file.write(f'Coordinates: {coordinates[0]:.6f} {coordinates[1]:.6f} {coordinates[2]:.6f}\n')

# Example usage
# write_input_data(handles)