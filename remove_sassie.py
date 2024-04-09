import os

from distutils.sysconfig import get_python_lib
lib_path = get_python_lib()
print (lib_path)

rmst = 'rm -Rf '+lib_path+'/sassie*'
result = os.popen(rmst).readlines()
for line in result: print (line)

rmst = 'rm -Rf ./build/'
result = os.popen(rmst).readlines()
for line in result: print (line)
      
