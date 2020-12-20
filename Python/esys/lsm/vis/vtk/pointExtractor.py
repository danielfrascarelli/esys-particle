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

import vtk as kwvtk
from esys.lsm.vis import core

class PointExtractor(core.PointExtractor):
    def __init__(
        self,
        pointMap = lambda dataRecord: dataRecord.getPoint()
    ):
        core.PointExtractor.__init__(self, pointMap)

    def getVtkPoints(self, data):
        vtkPoints = kwvtk.vtkPoints()
        for dataRecord in data:
            vtkPoints.InsertNextPoint(self.getPoint(dataRecord))
        return vtkPoints
