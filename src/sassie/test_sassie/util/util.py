import array

'''
set the abnormal float numbers for unittest
'''
import numpy

import os
floattype=os.environ['SASSIE_FLOATTYPE']


# Convert a number to a floattype ft
def num_to_floattype(a, ft):
    b = []
    b.append(a)
    if ft=='float32':
       result = array.array('f',b)[0]
    elif ft=='float64':
       result = array.array('d',b)[0]
    return result


# Recursively convert a list to a floattype type list
def list_to_floattype(lo, ft):
    try:
        l = list(lo)
    except:
        raise "ListConvertError"
    ltmp = []
    for i in range(len(l)):
        if isinstance(l[i], list):
           ltmp.append(list_to_floattype(l[i], ft))
        else:
           ltmp.append(num_to_floattype(l[i],ft))
    return ltmp


# Recursively format-print a list
import sys
def printfl(l, fmt='%.3f'):
   for i in l:
      if isinstance(i, (list,numpy.ndarray)):
         sys.stdout.write('[')
         printfl(i, fmt)
         sys.stdout.write('\b\b], ')
      else:
         sys.stdout.write(fmt%float(i)+', ')


NAN = float('nan')
INF = float('inf')
HUGE = float(numpy.finfo(floattype).max)
TINY = float(numpy.finfo(floattype).tiny)
ZERO = num_to_floattype(0.0, floattype)


if __name__=="__main__":
   ft = 'float32'
   a=[[1,2],[3,4,[5]],[[6,[7,8]],6],0]
   na=[]
   list_to_floattype(a,ft)
   print a,'\n',na

