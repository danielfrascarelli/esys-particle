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

"""
Constructive solid geometry classes: L{Union}.
"""

from .modifier import Modifiable

class Union(Modifiable):
    """
    Union of objects.
    """
    def __init__(self, *args):
        Modifiable.__init__(self)
        self.objList = []
        for arg in args:
            self.append(arg)

    def __iter__(self):
        return iter(self.objList)

    def __getattr__(self, attr):
        return getattr(self.objList, attr)

    def writeBegin(self, f):
        f.write("\nunion {\n")

    def writeEnd(self,f):
        f.write("\n}")

    def writeObjects(self,f):
        for obj in iter(self):
            obj.writeSdl(f)
            f.write("\n")

    def writeSdl(self, f):
        self.writeBegin(f)
        self.writeObjects(f)
        Modifiable.writeSdl(self, f)
        self.writeEnd(f)

