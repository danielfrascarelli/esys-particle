0#############################################################
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

"""
Defines various classes which control the parameters (arguments)
for running of a povray process.
"""
import string
import os
import tempfile

import esys.lsm.Logging
from   esys.lsm.vis.core.exception import raiseNotImplemented
from   esys.lsm.vis                import core
from   esys.lsm.util               import process

logger = esys.lsm.Logging.getLogger("esys.lsm.vis.povray")

class Option:
    def __init__(self, name):
        self.name  = name

    def getName(self):
        return self.name

    def getArgList(self):
        raiseNotImplemented("Implement this method in derived class.")

class SingleOption(Option):
    def __init__(self, name, arg, val = None):
        Option.__init__(self, name)
        self.arg   = arg
        self.value = val

    def getArg(self):
        return self.arg

    def getValue(self):
        return self.value

    def setValue(self, val):
        self.value = val

    def getArgList(self):
        if (self.getValue() != None):
            return ["+" + self.getArg() + str(self.getValue())]
        return []

class BoolOption(SingleOption):
    def __init__(self, name, arg, val = None):
        SingleOption.__init__(self, name, arg, val = None)

    def getValue(self):
        return bool(SingleOption.getValue(self))
                                            
    def getArgList(self):
        if (self.getValue() != None):
            argChar = "-"
            if (self.getValue()):
                argChar = "+"
            return [argChar + self.getArg()]

        return []

class BoolValOption(SingleOption):
    def __init__(self, name, arg, val = None):
        SingleOption.__init__(self, name, arg, val)

    def getArgList(self):
        if (self.getValue() != None):
            if (len(str(self.getValue())) > 0):
                return ["+" + self.getArg() + str(self.getValue())]
            else:
                return ["-" + self.getArg()]

        return []

class ListOption(SingleOption):
    def __init__(self, name, arg, val = None):
        SingleOption.__init__(self, name, arg)
        self.setValue([])

    def appendValue(self, val):
        self.getValue().append(val)

    def getArgList(self):
        return \
            [
                "+" + self.getArg() + str(val) \
                for val in self.getValue()
            ]

class DisplayPauseOption(Option):
    def __init__(self, name):
        Option.__init__(self, name)
        self.display = BoolOption("Display", "D", None)
        self.pause   = BoolOption("Pause", "P", None)
        self.noPause = BoolOption("Pause", "P", False)

    def setOffScreen(self, offScreen):
        self.display.setValue(not offScreen)

    def getOffScreen(self):
        return (not self.display.getValue())

    def setPause(self, pause):
        self.pause.setValue(pause)

    def getPause(self):
        return self.pause.getValue()

    def getArgList(self):
        argList = self.display.getArgList()
        if (self.getOffScreen()):
            argList += self.noPause.getArgList()
        else:
            argList += self.pause.getArgList()
        return argList

_formatCharDict = dict()
_formatCharDict[str(core.PNG)] = "N"
_formatCharDict[str(core.PNM)] = "P"

from esys.lsm.util import InstallInfo

class PovProcess(process.Process):
    def __init__(self, argList, combineStdOutAndStdErr=True):
        process.Process.__init__(
            self,
            InstallInfo.getPovrayExePath(),
            argList,
            combineStdOutAndStdErr,
            printOutput=False
        )

class PovSdlInput:
    def __init__(self):
        self.tempFile = tempfile.mkstemp()
        """
        Originally this was opened for writing as a binary file ('wb').
        But in object.py writes (e.g. class Background) are as strings.
        Under Python 2.x this appears to have worked because the strings
        were implicitly interpreted as ascii encoded.  Under Python 3.x
        unicode is the default encoding and so text writes to a binary
        file work not.
        """
        self.file = os.fdopen(self.tempFile[0], 'w')
        self.fileName = self.tempFile[1]

    def __del__(self):
        self.close()
        os.remove(self.getFileName())

    def getFile(self):
        return self.file

    def getFileName(self):
        return self.fileName

    def close(self):
        if (self.file != None):
          self.file.close()
          self.file = None

class PovRenderer:
    def __init__(self):
        self.nameOptDict = dict()
        self.initialiseOptionDict()

    def initialiseOptionDict(self):
        self.nameOptDict["Quality"]          = SingleOption("Quality", "Q")
        self.nameOptDict["Quality"].setValue(9)
        self.nameOptDict["Antialias"]        = SingleOption("Antialias", "A")
        self.nameOptDict["Antialias"].setValue(0.1)
        self.nameOptDict["DisplayPause"]     = DisplayPauseOption("DisplayPause")
        self.nameOptDict["Height"]           = SingleOption("Height", "H")
        self.nameOptDict["Width"]            = SingleOption("Width", "W")
        self.nameOptDict["Input_File_Name"]  = SingleOption("Input_File_Name", "I")
        self.nameOptDict["Output_To_File"]   = \
            BoolValOption("Output_To_File", "F")
        self.nameOptDict["Output_File_Name"] = SingleOption("Output_File_Name", "O")
        self.nameOptDict["Library_Path"]     = ListOption("Library_Path", "L")
        self.nameOptDict["All_Console"]      = \
            BoolValOption("All_Console", "G", "A")

    def getOption(self, optionName):
        return self.nameOptDict[optionName]

    def setOffScreen(self, doOffScreen):
        self.getOption("DisplayPause").setOffScreen(doOffScreen)

    def getOffScreen(self):
        return self.getOption("DisplayPause").getOffScreen()

    def setInteractive(self, interactive):
        self.getOption("DisplayPause").setPause(interactive)

    def getInteractive(self):
        return self.getOption("DisplayPause").getPause()

    def getImageFormat(self, fileName):
        return core.getFormatFromExtension(fileName)

    def getImageFormatChar(self, imageFormat):
        return _formatCharDict[str(imageFormat)]

    def setFileName(self, fileName, imageFormat=None):
        self.getOption("Output_File_Name").setValue(fileName)
        if (fileName != None):
            if (imageFormat == None):
                imageFormat = self.getImageFormat(fileName)
            self.getOption("Output_To_File").setValue(
                self.getImageFormatChar(imageFormat)
            )
        else:
            self.getOption("Output_To_File").setValue("")

    def getFileName(self):
        return self.getOption("Output_File_Name").getValue()

    def setSize(self, size):
        self.getOption("Width").setValue(size[0])
        self.getOption("Height").setValue(size[1])

    def getSize(self):
        return \
          [
              self.getOption("Width").getValue(),
              self.getOption("Height").getValue()
          ]

    def clear(self):
        pass

    def getArgList(self):
        argList = []
        for opt in list(self.nameOptDict.values()):
            argList += opt.getArgList()
        return argList

    def openInput(self):
        self.sdlInput = PovSdlInput()
        return self.sdlInput.getFile()

    def getInput(self):
        return self.sdlInput.getFile()

    def closeInput(self):
        self.sdlInput.close()
        self.getOption("Input_File_Name").setValue(self.sdlInput.getFileName())
        povproc = PovProcess(self.getArgList())
        logger.info("Running process: " + povproc.getCommandLine())
        povproc.run()
        logger.info("povray output = " + str(povproc.getStdOutLineList()))

