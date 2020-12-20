"""
Provides install path infomation.
"""
import os
from esys.lsm.util.pathSearcher import PathSearcher

installDir   = "/home/daniel/Documents/fing/esys-particle/src/danielfrascarelli-git/esys-particle"
binDir       = os.path.join(installDir, "bin")
libDir       = os.path.join(installDir, "lib")

pythonPkgDir = "/home/daniel/Documents/fing/esys-particle/src/danielfrascarelli-git/esys-particle/lib/python2.7/site-packages"
esysPkgDir   = os.path.join(pythonPkgDir, "esys")
lsmPkgDir    = os.path.join(esysPkgDir, "lsm")

pkgName      = "ESyS-Particle"
version      = "2.3.5"
pkgHomePageUrl = "https://launchpad.net/esys-particle/"
pkgDataDir   = "/home/daniel/Documents/fing/esys-particle/src/danielfrascarelli-git/esys-particle/share/esys-particle"
povrayExe    = "no"

_haveVtk     = False
_havePovray  = False

def getPovrayExePath():
    """
    Attempts to return the absolute path of the "povray" executable
    using the "PATH" environment variable. If the exe can't be found
    on the "PATH" then this function returns the "povray" path which
    was found during installation. This function is a workaround for
    for the SGI MPT mpirun, which seems to alter the user "PATH"
    environment.
    """
    absPath=PathSearcher().find("povray")
    if ((absPath == None) or (absPath == "")):
        absPath = povrayExe

    return absPath

def getDataFilePath(dataFileName):
    """
    Returns path for specified data file. Looks on path
    C{L{pkgDataDir}:Data:.}
    """
    return PathSearcher(pkgDataDir+":Data:.").which(dataFileName)

def haveVtk():
    return _haveVtk

def havePovray():
    return _havePovray

