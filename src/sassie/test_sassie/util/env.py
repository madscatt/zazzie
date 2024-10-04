import os
import math



def determine_float_type():
   '''
      determine whether it will be a 32/64-bit number for the float type
   '''
   a=1.e40
   if math.isinf(a):
      return 'float32'
   # set up a pseudo-32bit machine for testing
   elif a==1.e40:
      return 'float64'
   # Need a reasonable default
   else:
      return 'float'


if os.environ.has_key('SASSIE_LARGETEST'):
   print '\nKEY CONFLICT IN os.envrion of SASSIE_LARGETEST!\nWILL QUIT!'
   exit()
else:
   os.environ['SASSIE_LARGETEST']='n'


if os.environ.has_key('SASSIE_HUGETEST'):
   print '\nKEY CONFLICT IN os.envrion of SASSIE_HUGETEST!\nWILL QUIT!'
   exit()
else:
   os.environ['SASSIE_HUGETEST']='n'


# Genrate huge dcd files if SASSIE_HUGETEST is specified
# IMPORTANT, SASSIE_HUGETEST only implies the huge rna dcd files generated below will be used
if os.environ['SASSIE_HUGETEST']=='y':
	from sassie.core_testing.util import generate_huge_dcd_onthefly
	generate_huge_dcd_onthefly.generate_huge_dcd()


if os.environ.has_key('SASSIE_FLOATTYPE'):
   print '\nKEY CONFLICT IN os.envrion of SASSIE_FLOATTYPE!\nWILL QUIT!'
   exit()
else:
   os.environ['SASSIE_FLOATTYPE']= determine_float_type()



"""
while ft!='' and ft!='1' and ft!='2' and ft!='3':
   ft = raw_input('\nPlease choose the floattype for your test:\n',\
      '1: automatic\n2: float32\n3: float64\n'\
      'Please select-->(default automatic')

if ft=='' or ft=='1':
   floattype = determine_float_type()
elif ft=='2':
   floattype = 'float32'
elif ft=='3':
   floattype = 'float64'
"""

"""
floattype = 'float'
ft = [item for item in sys.argv if 'floattype' in item]
if len(ft)>0:
   if len(ft)>1:
      print 'FLOATTYPE SPECIFIED MORE THAN ONCE.\nWILL ONLY USE THE FIRST PLACE OF SPECIFICATION!\n'
   floattype = ft[0].split('=')[1]
else:
   floattype = determine_float_type()


largetest = False
if '-largetest' in sys.argv:
   largetest = True
else:
   largetest = False


hugetest = False
if '-hugetest' in sys.argv:
   hugetest = True
else:
   hugetest = False
"""
