import numpy as np
import math
import random

import warnings

warnings.filterwarnings('ignore')

def rotation_matrix(A,B):
# a and b are in the form of np array

   ax = A[0]
   ay = A[1]
   az = A[2]

   bx = B[0]
   by = B[1]
   bz = B[2]

   au = A/(np.sqrt(ax*ax + ay*ay + az*az))
   bu = B/(np.sqrt(bx*bx + by*by + bz*bz))

   R=np.array([[bu[0]*au[0], bu[0]*au[1], bu[0]*au[2]], [bu[1]*au[0], bu[1]*au[1], bu[1]*au[2]], [bu[2]*au[0], bu[2]*au[1], bu[2]*au[2]] ])


   return(R)

attempts = 0
failures = 0
for attempt in xrange(10000000):
    a = np.array([1,0,0])
    b = np.array([0,1,0])
    b = np.array([random.randint(0,10),random.randint(0,10),random.randint(0,10)])

    try:
        rot = rotation_matrix(a,b)
    except:
        break

    c = rot.dot(a)
    try:
        cu = c/np.linalg.norm(c)
        bu = b/np.linalg.norm(b)
    except:
        break

    attempts += 1
    dot_product = cu.dot(bu)

    if dot_product == 1.0:
	    pass
	    #print 'GOOD'
	    #print 'cu == bu'
	    #print 'cu = ', cu
	    #prink 'bu = ', bu

    elif np.abs(np.sum(bu - cu)) > 0.001:
	    print 'UH OH!'
	    print 'cu != bu'
	    print 'bu - cu = ', bu - cu
	    failures += 1

print 'attempts = ',attempts
print 'failures = ',failures
