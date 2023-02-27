import numpy as np
import dna_overlap
npart = 9000
coor = np.zeros((npart,3))
for i in range(npart):
    coor[i]=np.random.rand(3)
distances = dna_overlap.distances(coor,np.zeros((npart,npart)))
print(coor)
print(distances)
