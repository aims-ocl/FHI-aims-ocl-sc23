'''
Created on Oct 30, 2013

@author: newhouse
'''
from profilestats import profile
import time

from core.database import delete_db, initialize_tables
from core.run_GA import RunGA
import core.run_GA

def profile():
    delete_db()
    initialize_tables()
    
    ga = RunGA()
    ga.scan_initial_pool()
#     ga.start()

profile()