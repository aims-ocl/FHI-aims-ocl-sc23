'''
Created on Aug 5, 2013

@author: newhouse
'''

# adding src directory to path at runtime
import os
import shutil
import subprocess
import sys
import time
# add source directory to python path
src_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.append(src_dir)

from core import user_input
from core.file_handler import *
from core.kill import set_kill




def main(*args, **kwargs):
    # read command line arguments
    argv = args[0][1:]
    
    if 'clean' in argv or 'reset' in argv:
        clean()
        return
    elif 'datatools' in argv:
        from core import data_tools
        data_tools.main(argv)
        return
    elif 'kill' in argv:
        set_kill()
        return
    else:
        run_GA(argv)    
    
    
def fill_initial_pool():
    set_progress('starting')

    os.system('python ' + os.path.join(src_dir, 'core', 'fill_initial_pool.py'))
    # replace with following for thnec
    # os.system('qsub fill_initial_pool_jobs')
    
def run_GA(argv):
    minor_version = sys.version_info[1]
    ui = user_input.get_config()
    number_of_replicas = ui.get_eval('run_settings', 'number_of_replicas')
    try: stoic_filename = argv[0]
    except: stoic_filename = '' 
    for i in range(number_of_replicas):
        if minor_version == 7: os.system('python ' + os.path.join(src_dir, 'core', 'run_GA.py ') + stoic_filename)  # laptop
        if minor_version == 6: os.system('qsub jobs')  # thnec        
    
def clean():
    print 'resetting environment'
    directory_to_remove = [tmp_dir, structure_dir]
    files_to_remove = [output_file]
    p = subprocess.Popen(['rm *.out'], cwd=cwd, shell=True)
    p.wait()
    # tmp index is to keep track of replica number
    for directory in directory_to_remove:
        if os.path.exists(directory): rmtree(directory)
    for rmfile in files_to_remove:
        if os.path.exists(rmfile): os.remove(rmfile)
    return


if __name__ == '__main__':
    main(sys.argv)
