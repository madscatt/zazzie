import filecmp


def cmp_skip(f1, f2, skip=[]):
   '''
   compare two text files by skipping some lines
   '''
   print 'skipping ',skip
   with open(f1,'r') as fp1, open(f2,'r') as fp2:
      count = 0
      while True:
         count += 1
         line1 = fp1.readline()
         line2 = fp2.readline()
         if (count in skip):
            continue
         elif (not line1):
            return True
         elif (line1 != line2):
            print '\ndiffer beginning at:\n',line1, '\n', line2, '\n'
            return False
      return True

def cmp_dcd_skip_date(f1,f2,start=179,end=202):
   '''
   compare two dcd files
   '''
   fp1 = open(f1, 'rb')
   fp2 = open(f2, 'rb')
   fcont1 = fp1.read()
   fcont2 = fp2.read()
   l = len(fcont1)
   if l!=len(fcont2):
      print "File length does not match"
      return False
   else:
      for i in range(l):
         #print fcont1[i],fcont2[i]
         if i>start and i<end:
            continue
         else:
            if fcont1[i]!=fcont2[i]:
               return False
   return True

if __name__ == '__main__':
   cmp_skip('a','b',[4, 59])
   print cmp_dcd_skip_date('c.dcd','d.dcd')

