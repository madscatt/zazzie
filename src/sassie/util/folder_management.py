import os,glob,locale,re,shutil

def back_up_existing_module_folder(runname,module):
    folders = glob.glob(os.path.join(runname,module+'*'))
    if os.path.join(runname,module) not in folders:
        return
    max_num = 0
    for folder in folders:
        folder_name = os.path.basename(folder)
        if re.compile('^'+module+'_\d+$').match(folder_name):
            num = locale.atoi(folder_name[len(module)+1:])
            if num>max_num:
                max_num = num
    shutil.move(os.path.join(runname,module), os.path.join(runname,module+'_%d'%(max_num+1)))

    return

if __name__ == '__main__':

    path = './'
    runname = 'run_fred'
    module = 'happy_feet'
    back_up_existing_module_folder(runname,module)


