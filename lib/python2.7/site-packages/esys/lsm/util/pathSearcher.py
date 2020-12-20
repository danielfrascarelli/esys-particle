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

import os

class PathSearcher:
    
    def __init__(self, pathString=os.environ["PATH"], pathDelimiter=":"):
        self.searchList = str.split(pathString, pathDelimiter)
        
    def contains(self, fileName):
        for dir in self.searchList:
            if os.path.exists(os.path.join(dir, fileName)):
                return True
        return False

    def which(self, fileName):
        for dir in self.searchList:
            fileNamePath = os.path.join(dir, fileName)
            if os.path.exists(fileNamePath):
                return fileNamePath
        return fileName
        
    def find(self, fileName):
        for dir in self.searchList:
            fileNamePath = os.path.join(dir, fileName)
            if os.path.exists(fileNamePath):
                return fileNamePath
        return ""
