#!/usr/bin/python3
"""
This module is part of the revised regression tests for FHI-aims and provides:
    -AbstractRemote: abstract class for connectors to remote servers, prepare
                     the necessary information for the queueing system on the
                     destination server and commit your job there

Note that this is a backend module. If came here to add a new regression test,
you are searching the wrong place. Please read the documentation in the main
folder for guidelines.
"""

import tarfile
import os.path
import subprocess

class AbstractRemote(object): #pylint: disable = R0921
    """this is the base class for any remote server connector"""

    #--------------------------------------------------------------
    #constructor function(s)
    #--------------------------------------------------------------

    def __init__(self, data, localworkdir):
        self._sshtarget = data["ssh"]
        self._workdir = data["remoteworkdir"]
        self._archivename = data["archivename"]
        self._data = self._parse_config(data["remoteconfig"])
        self._locworkdir = localworkdir
        self._globdata = data
        self._stream = None

    #--------------------------------------------------------------
    #protocol functions
    #--------------------------------------------------------------

    def __enter__(self):
        path = os.path.join(self._locworkdir, self._data["jobscriptname"])
        self._stream = open(path, "w")
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._stream.close()

    #--------------------------------------------------------------
    #properties
    #--------------------------------------------------------------

    #--------------------------------------------------------------
    #protected functions
    #--------------------------------------------------------------

    def _connect(self):
        """copy the archive and connect to the remote, using SSH"""
        path = os.path.join(self._locworkdir, self._archivename+".tar.gz")
        remotepath = os.path.join(self._workdir, self._archivename+".tar.gz")
        print("ssh-connecting to create remote working directory...")
        args = ["ssh", self._sshtarget, 'mkdir -p "%s"'%(self._workdir)]
        retval = subprocess.call(args)
        if retval != 0:
            msg = "Failed to create remote working directory!"
            raise subprocess.CalledProcessError(retval, args, msg)
        print("transferring archive to remote server...")
        args = ["scp", path, "%s:%s"%(self._sshtarget, remotepath)]
        retval = subprocess.call(args)
        if retval != 0:
            msg = "Archive could not be transferred to remote server"
            raise subprocess.CalledProcessError(retval, args, msg)
        print("ssh-connecting to remote server to queue job...")
        args = ["ssh", self._sshtarget]
        path = os.path.join(self._locworkdir, "submit.sh")
        subprocess.call(args, stdin=open(path, "r"))

    def _generate_shell(self):
        """
        create the shell script with the necessary commands to queue the job on
        the remote server, the file should be named "submit.sh" to run with the
        default implementation of the _connect() function. Do not forget to
        close the (SSH-)connection at the end of the script.
        """
        raise NotImplementedError

    def _pack_archive(self):
        """
        pack the archive for transmission to the remote server, defaults to
        creating a tar.gz archive, but can be subclassed if another format is
        necessary"""
        archive = os.path.join(self._locworkdir, self._archivename+".tar.gz")
        with tarfile.open(archive, "w:gz") as tarball:
            tarball.add(self._locworkdir, self._archivename)

    def _parse_config(self, filename):
        """
        parse the given config file and store all relevant data in a dictionary
        which should be returned at the end of the function (to be stored on
        this object by the constructor). Note that a "jobscriptname" entry is
        required.
        """
        raise NotImplementedError

    def _write(self, string, newline=True):
        """write a line to the jobscript outstream while the file is open"""
        if newline:
            self._stream.write(string+"\n")
        else:
            self._stream.write(string)

    #--------------------------------------------------------------
    #public functions
    #--------------------------------------------------------------

    def finalize_preparation(self):
        """
        the final steps to do before the regression test is ready for submission
        """
        self._write("tar --gzip -cf ../%s_results.tar.gz reference test"%(
            self._archivename))

    @classmethod
    def typename(cls):
        """
        a tuple with the shorthand notation for this type and a short
        description
        """
        raise NotImplementedError

    def process_testcase(self, case):
        """make the necessary preparations for a TestCase object"""
        raise NotImplementedError

    def process_testsuite(self, suite):
        """make the necessary preparations for a TestSuite object"""
        raise NotImplementedError

    def process_regression(self, suite):
        """make the top-level preparations for the whole regression test"""
        raise NotImplementedError

    def submit(self):
        """submit everything relevant to the server"""
        self._pack_archive()
        self._generate_shell()
        self._connect()
