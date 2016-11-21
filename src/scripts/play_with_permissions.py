import os
''' make a directory '''
os.system('mkdir empty_folder')
''' see if you can read the directory '''
print os.access('empty_folder', os.R_OK)
''' make the directory un-readable'''
os.system('chmod a-r empty_folder')
''' see if you can read the directory '''
print os.access('empty_folder', os.R_OK)

''' make the directory readable'''
os.system('chmod a+r empty_folder')
''' remove the directory '''
os.system('rm -Rf empty_folder')
