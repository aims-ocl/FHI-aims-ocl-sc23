#!/usr/bin/python3
"""
This module is part of the revised regression tests for FHI-aims and provides:
    -LoadLevelerRemote - communication with remote cluster using IBM LoadLeveler

Note that this is a backend module. If came here to add a new regression test,
you are searching the wrong place. Please read the documentation in the main
folder for guidelines.
"""

import os.path

from .abstractremote import AbstractRemote

class LoadLevelerRemote(AbstractRemote):
    """
    this implementation of the AbstractRemote class is designed for clusters
    that run the Solaris Grid Engine as queueing system
    """

    #--------------------------------------------------------------
    #constructor function(s)
    #--------------------------------------------------------------

    #--------------------------------------------------------------
    #protocol functions
    #--------------------------------------------------------------

    #--------------------------------------------------------------
    #properties
    #--------------------------------------------------------------

    #--------------------------------------------------------------
    #protected functions
    #--------------------------------------------------------------

    def _generate_shell(self):
        with open(os.path.join(self._locworkdir, "submit.sh"), "w") as shell:
            shell.write("#!/bin/bash\n")
            shell.write("cd %s\n"%self._workdir)
            shell.write("tar -xf %s.tar.gz\n"%self._archivename)
            shell.write("cd %s\n"%self._archivename)
            shell.write("llsubmit regression.sge\n")
            shell.write("exit\n")

    def _parse_config(self, filename):
        data = {"jobscriptname":"regression.sge", "template":filename}
        return data

    #--------------------------------------------------------------
    #public functions
    #--------------------------------------------------------------

    @classmethod
    def typename(cls):
        return "IBM", "IBM LoadLeveler"

    def process_regression(self, suite):
        with open(self._data["template"], "r") as instream:
            for line in instream:
                self._write(line, newline=False)
        if not os.path.isdir(suite.reference):
            self._write('aims_ref="%s"'%(suite.reference))
        self._write('aims_test="%s"'%(suite.testsubject))
        basefolder = os.path.join(self._workdir, self._archivename)
        self._write('basefolder="%s"\n'%basefolder)
        self._write('cd $basefolder\n')

    def process_testcase(self, case):
        if not os.path.isdir(case.reference):
            basepath = os.path.join("reference", case.folder_case)
            for folder in case.subfolders:
                path = os.path.normpath(os.path.join(basepath, folder))
                self._write('cd "%s"'%path)
                self._write('%s $aims_ref 1> aims.out 2> aims.err'%(
                    " ".join(self._globdata["command"])))
                self._write('cd $basefolder')
        basepath = os.path.join("test", case.folder_case)
        for folder in case.subfolders:
            path = os.path.normpath(os.path.join(basepath, folder))
            self._write('cd "%s"'%path)
            self._write('%s $aims_test 1> aims.out 2> aims.err'%(
                " ".join(self._globdata["command"])))
            self._write('cd $basefolder')
        self._write("")

    def process_testsuite(self, suite):
        pass
