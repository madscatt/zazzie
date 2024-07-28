import numpy as np
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
                value = handles[key]
                if isinstance(value, (list, np.ndarray)):
                    file.write(f"{key}:\n")
                    for item in value:
                        if isinstance(item, (list, np.ndarray)):
                            formatted_item = " ".join(
                                str(int(x) + 1) if isinstance(x, int) else str(x) for x in item
                            )
                            file.write(formatted_item + "\n")
                        else:
                            file.write(str(int(item) + 1) if isinstance(item, int) else str(item) + "\n")
                else:
                    file.write(f"{key}: {pprint.pformat(value)}\n")
            else:
                file.write(f"{key}: Key not found in handles\n")

    return

print('DONE writing input data for calcI_SlblAtoms')

# Example usage
# handles = {
#     'nEnCalc': 10,
#     'nDCalc': 5,
#     'nSCalc': 3,
#     'nlblTypes': 2,
#     'nLabels': 7,
#     'labelData': [[1, 2, 3], [4, 5, 6]],
#     'atomData': np.array([[1, 2, 3], [4, 5, 6]]),
#     'deltaD': 0.1,
#     'fPairCalc': [0.2, 0.3],
#     'actDMax': 1.0
# }
# write_input_data(handles)