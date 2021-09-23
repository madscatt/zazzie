

#input_list =  1,2; 3,4; 5,6
#input_list.split =  ['1,2', ' 3,4', ' 5,6']
#nested_array =  [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]

import string
import locale

start_list = '1.1,2.0; 3.2,4.1; -5.3,6.2'
print 'start_list = ', start_list

svariables={}

svariables['delta_rho'] = (start_list,'nested_float_array')

for key in svariables:
    if (svariables[key][1] == 'nested_float_array'):
        input_list = svariables.get(key)[0]
        print 'input_list: ', input_list
#        print 'value[0]: ', value[0]

    lin = input_list.split(';')

    print 'input_list.split = ', lin

    nested_array = []
    for item in lin:
        new_item = item.replace(' ', '')
        print 'new item: ', new_item
        new_item = new_item.split(',')
        print 'new item 2: ', new_item
        temp = []
        for value in new_item:
            float_item = locale.atof(value)
            temp.append(float_item)

        nested_array.append(temp)

    print 'nested_array = ', nested_array
    print 'len(nested_array): ', len(nested_array)
    print 'type(nested_array): ', type(nested_array)