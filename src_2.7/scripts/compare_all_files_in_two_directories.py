from __future__ import print_function
import filecmp
import os

"""
    Simple file comparison check between all files in two directories.
    This script does not search recursively in the provided directories.

    JEC: 07/22/16

"""

directory_1 = '/Users/curtisj/git_working_copies/zazzie/src/sassie/calculate/sascalc/sascalc_library/cpp_and_cuda_buildBlock/src'
directory_2 = '/Users/curtisj/git_working_copies/zazzie/src/core_libraries_1.0/cpp_and_cuda/sascalc_library/cpp_and_cuda/src'

directory_1 = '/Users/curtisj/git_working_copies/zazzie/src/sassie/calculate/sascalc/sascalc_library/cpp_extension'
directory_2 = '/Users/curtisj/git_working_copies/zazzie/src/core_libraries_1.0/cpp_and_cuda/sascalc_library/cpp_extension'

# Determine the items that exist in both directories
d1_contents = set(os.listdir(directory_1))
d2_contents = set(os.listdir(directory_2))
common = list(d1_contents & d2_contents)
common_files = [ f 
                for f in common 
                  if os.path.isfile(os.path.join(directory_1, f))
                  ]

print('Common files:', common_files)

# Compare the directories
match, mismatch, errors = filecmp.cmpfiles(directory_1, directory_2, common_files)

#.report_full_closure()

print()
print('Match:', match)
print('Mismatch:', mismatch)
print('Errors:', errors)

[print() for x in xrange(3)]

filecmp.dircmp(directory_1, directory_2).report_full_closure()

print()
