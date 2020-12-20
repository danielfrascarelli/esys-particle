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
import os.path

def removeDirs(dir):
    """
        Recursively removes directory and all contents and
        subdirectories.
    """
    for f in os.listdir(dir):
        absFile = os.path.join(dir, f)
        if (os.path.isdir(absFile)):
            removeDirs(absFile)
        else:
            os.remove(absFile)
    os.rmdir(dir)
