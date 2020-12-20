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

class SurfaceData(core.SurfaceData):
    def __init__(self, data, pointExtractor):
        core.SurfaceData.__init__(self, data, pointExtractor)

    def getVtkPolyData(self):
        vtkPolyData = kwvtk.vtkPolyData()
        vtkPolyData.SetPoints(
            self.getPointExtractor().getVtkPoints(self.getData())
        )
        return vtkPolyData

    def getVtkDelaunay(self):
        vtkDelny = kwvtk.vtkDelaunay2D()
        vtkDelny.SetInput(self.getVtkPolyData())
        return vtkDelny

    def getVtkPolyDataNormals(self):
        vtkNormals = kwvtk.vtkPolyDataNormals()
        vtkNormals.SetInput(self.getVtkDelaunay().GetOutput())
        return vtkNormals

    def getVtkMapper(self):
        vtkMapper = kwvtk.vtkPolyDataMapper()
        vtkMapper.SetInput(self.getVtkPolyDataNormals().GetOutput())
        return vtkMapper

    def getActor(self):
        vtkActor = kwvtk.vtkActor()
        vtkActor.SetMapper(self.getVtkMapper())

        return vtkActor

    def addActor(self, addTo):
        addTo.AddActor(self.getActor())


