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

import logging

strLevelDict = dict()
strLevelDict["INFO"]     = logging.INFO
strLevelDict["WARNING"]  = logging.WARNING
strLevelDict["ERROR"]    = logging.ERROR
strLevelDict["DEBUG"]    = logging.DEBUG
strLevelDict["CRITICAL"] = logging.CRITICAL

def getLevel(levelString):
    return strLevelDict[str.upper(levelString)]

def getStringLevelList():
    return list(strLevelDict.keys())

getLogger = logging.getLogger
basicConfig = logging.basicConfig
