import os

import numpy as np
import sassie.sasmol.sasmol as s
import sassie.sasmol.dcdio as dcdio


def generate(fin='rna.pdb', fout='rna.dcd', frames=1000):
  path = os.path.dirname(os.path.realpath(__file__))
  fin = os.path.join(path,'../','data','pdb_common',fin)
  if not os.path.exists(fin):
    raise Exception,fin+' does not exist while generating huge rna dcd files'
    return
  #fout = os.path.join(path,'../','data','dcd_common',fout)
  fout = os.path.join('/tmp/',fout)
  if os.path.exists(fout):
    return
  o=s.SasMol(0)
  o.read_pdb(fin)
  outfile = dcdio.open_dcd_write(fout)
  nset = frames
  natoms = o._coor[0,:,0].shape[0]
  istart = 0 ; nsavc = 1 ; delta = 1.0	
  o.write_dcd_header(outfile,frames)
  tx=o.coor()[0,:,0].astype(np.float32)	
  ty=o.coor()[0,:,1].astype(np.float32)	
  tz=o.coor()[0,:,2].astype(np.float32)	
  for i in range(frames):
    dcdio.write_dcdstep(outfile,tx,ty,tz,i+1)
  o.close_dcd_write(outfile)

def generate_huge_dcd():
  fin = 'rna.pdb'
  data = {'rna-0.8g.dcd':6250, 'rna-1.0g.dcd':7813, 'rna-1.2g.dcd':9375, 'rna-2.0g.dcd':15625, 'rna-3.2g.dcd':25000, 'rna-6.4g.dcd':50000}
  for (value,key) in zip(data.keys(),data.values()):
    generate(fin,value,key)

if __name__=='__main__':
  import datetime
  start = datetime.datetime.now()
  generate_huge_dcd()
  end = datetime.datetime.now()
  print 'time used (seconds): ',(end-start).seconds
