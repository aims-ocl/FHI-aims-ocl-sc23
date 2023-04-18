'''
Created on Nov 29, 2013

@author: newhouse
'''

from file_handler import *

def set_unkill(): write_data(tmp_dir, 'kill.dat', 'live')
def set_kill(): write_data(tmp_dir, 'kill.dat', 'kill')
def get_kill(): return read_data(tmp_dir, 'kill.dat')

if __name__ == '__main__':
    set_kill()