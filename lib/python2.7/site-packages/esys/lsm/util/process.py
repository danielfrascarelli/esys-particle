#############################################################
##                                                         ##
## Copyright (c) 2003-2017 by The University of Queensland ##
## Centre for Geoscience Computing                         ##
## http://earth.uq.edu.au/centre-geoscience-computing      ##
##                                                         ##
## Primary Business: Brisbane, Queensland, Australia       ##
## Licensed under the Open Software License version 3.0    ##
## http://www.apache.org/licenses/LICENSE-2.0              ##
##                                                         ##
#############################################################

from esys.lsm.util.pathSearcher import *

import esys.lsm.Logging

import sys
import subprocess
import threading

def getLogger():
    return esys.lsm.Logging.getLogger("esys.lsm.util")

class ThreadedWriter(threading.Thread):
    def __init__(self, readFrom=None, writeToList=[]):
        threading.Thread.__init__(self)
        self.readFrom = readFrom
        self.writeToList = writeToList
        self.lineList = []

    def run(self):
        if (self.readFrom != None):
          try:
            line = self.readFrom.readline()
            """
            Writing a binary line to a text stream worked under
            Python 2.x without explicit decode, probably because
            the default encoding was ascii.  But under Python 3.x
            with unicode an explicit decode is necessary.
            """
            if sys.version_info[0] > 2:
              line = bytes.decode(line)
          except:
            line = ""
          while (line != ""):
            for f in self.writeToList:
              f.write(line)
            self.lineList.append(line)
            try:
              line = self.readFrom.readline()
              if sys.version_info[0] > 2:
                line = bytes.decode(line)
            except:
              line = ""

    def getLineList(self):
        self.join()
        return self.lineList

class Process:
    def __init__(
        self,
        exePath,
        exeArgList=None,
        combineStdOutAndStdErr=True,
        printOutput=False
    ):
        self.exePath                = PathSearcher().which(exePath)
        if (exeArgList == None):
            exeArgList = []
        self.exeArgList             = exeArgList
        self.printOutput            = printOutput
        self.combineStdOutAndStdErr = combineStdOutAndStdErr
        self.popen                  = None
        self.stdOutWriter           = None
        self.stdErrWriter           = None

    def getExePath(self):
        return self.exePath

    def getExeArgList(self):
        return self.exeArgList

    def getCommandLine(self):
        return \
            str(self.getExePath()) + " " + \
            str.join(" ", (list(map(str, self.getExeArgList()))))

    def wait(self):
        try:
            getLogger().info(
                "Waiting on process with id=" + \
                str(self.popen.pid)   + \
                ", cmd="     + \
                self.getCommandLine()
            )
            if sys.version_info[0] > 2:
              return self.popen.communicate()
            else:
              return self.popen.wait()
        except:
            getLogger().error(
                "Exception while waiting for cmd=" +\
                self.getCommandLine()
            )
            raise

    def runNoWait(self):
        if (self.printOutput):
            stdOutList = [sys.stdout]
            stdErrList = [sys.stderr]
        else:
            stdOutList = []
            stdErrList = []

        try:
            getLogger().info(
                "Creating process with cmd=" + self.getCommandLine()
            )
            if (self.combineStdOutAndStdErr):
                self.popen = \
                    subprocess.Popen(
                        self.getCommandLine().split(),
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        close_fds=True
                    )
            else:
                self.popen = \
                    subprocess.Popen(
                        self.getCommandLine().split(),
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        close_fds=True
                    
                    )
            self.stdOutWriter = ThreadedWriter(self.popen.stdout, stdOutList)
            self.stdErrWriter = ThreadedWriter(self.popen.stderr, stdErrList)
            getLogger().info("stdout = " + str(self.popen.stdout))
            getLogger().info("stderr = " + str(self.popen.stderr))
            getLogger().info("stdin  = " + str(self.popen.stdin))
            self.stdOutWriter.start()
            self.stdErrWriter.start()

        except:
            getLogger().error(
                "popen failed for cmd=" + self.getCommandLine() + "\n"
            )
            raise

    def run(self):
        self.runNoWait()
        return self.wait()

    def getStdInFile(self):
        return self.popen.stdin

    def getStdOutFile(self):
        return self.popen.stdout

    def getStdErrFile(self):
        return self.popen.stderr

    def isRunning(self):
        return (self.popen.poll() == -1)

    def getExitStatus(self):
        pollVal = self.popen.poll()
        if (pollVal == -1):
            return None
        return pollVal

    def getStdOutLineList(self):
        return self.stdOutWriter.getLineList()

    def getStdOutString(self):
        return str.join("", (self.stdOutWriter.getLineList()))

    def getStdErrLineList(self):
        return self.stdErrWriter.getLineList()

    def getStdErrString(self):
        return str.join("", (self.stdErrWriter.getLineList()))
