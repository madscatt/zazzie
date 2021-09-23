

#input_list =  1,2; 3,4; 5,6
#input_list.split =  ['1,2', ' 3,4', ' 5,6']
#nested_array =  [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]

import string
import locale

input_list = '1.1,2.0; 3.2,4.1; -5.3,6.2'
print 'input_list = ', input_list

input_list = input_list.split(';')

print 'input_list.split = ', input_list

nested_array = []
for item in input_list:
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