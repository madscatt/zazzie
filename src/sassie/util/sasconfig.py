import os

###     BEGIN SYSADMIN EDIT ###
###     BEGIN SYSADMIN EDIT ###
###     BEGIN SYSADMIN EDIT ###

#__arch__ = "cluster"
#__arch__ = "mac"
__arch__ = "linux"

__cuda__ = False
#__cuda__ = True

# total memory for single GPU card (4799 MiB) Tesla K20m
__gpu_memory__ = 4799

# __level__ = 'WARNING'
__level__ = 'DEBUG'



###     END SYSADMIN EDIT ###
###     END SYSADMIN EDIT ###
###     END SYSADMIN EDIT ###


if __arch__ == "cluster":

    installation_bin_path = ['share','apps','local','bin']
    __core_libraries_include__ = [os.path.join(os.path.sep, 'share', 'apps', 'local', 'core_libraries', 'include')]
    __core_libraries_lib__ = [os.path.join(os.path.sep, 'share', 'apps', 'local', 'core_libraries', 'lib')]

    if __cuda__:
        __cuda_path__ = os.path.join(os.path.sep, 'share', 'apps', 'local', 'cuda-6.5')
else:

    installation_bin_path = ['usr','local','bin']
    __core_libraries_include__ = [os.path.join(os.path.sep, 'usr', 'local', 'core_libraries', 'include')]
    __core_libraries_lib__ = [os.path.join(os.path.sep, 'usr', 'local', 'core_libraries', 'lib')]

    if __cuda__:
        __cuda_path__ = os.path.join(os.path.sep, 'usr', 'local', 'cuda-6.5')

__bin_path__ = os.path.join(os.path.sep,*installation_bin_path)



