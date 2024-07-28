import pprint

def write_input_data(handles, filename="input_data_for_calcI_Slblatoms.txt"):
    required_keys = [
        'nEnCalc', 'nDCalc', 'nSCalc', 'nlblTypes', 'nLabels', 'labelData', 
        'atomData', 'deltaD', 'fPairCalc', 'actDMax'
    ]
    
    with open(filename, 'w') as file:
        file.write("Input data for calcI_SlblAtoms:\n")
        for key in required_keys:
            if key in handles:
                file.write(f"{key}: {pprint.pformat(handles[key])}\n")
            else:
                file.write(f"{key}: Key not found in handles\n")

    return

print('DONE writing input data for calcI_SlblAtoms')

# Write input data to file before calculating label-atom scattering I(S)
#write_input_data_for_calcI_SlblAtoms(handles)

# Calculate the label-atom scattering I(S)
#handles, errorlg = calcI_SlblAtoms(handles)
#if errorlg:
#    raise Exception('Label-scatter intensity calc error')

"""
# Calculate the label-label scattering I(S)
"""