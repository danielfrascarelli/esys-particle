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

from __future__   import division
from esys.lsm.vis import core

import vtk as kwvtk

class EdgeData(core.EdgeData):
    def __init__(self, data, edgeExtractor):
        core.EdgeData.__init__(self, data, edgeExtractor)

    def getVtkData(self):
        return self.getEdgeExtractor().getVtkData(self.getData())

    def getVtkTube(self):
        vtkData = self.getVtkData()
        vtkTube = kwvtk.vtkTubeFilter()
        vtkTube.SetVaryRadiusToVaryRadiusByScalar()
        vtkTube.SetNumberOfSides(16)
        vtkTube.SetInput(vtkData)
        radiusRange = \
            vtkData.GetPointData().GetArray(
                self.getEdgeExtractor().getRadiusFieldName()
            ).GetRange(
                self.getEdgeExtractor().getRadiusFieldComponentIndex()
            )
        vtkTube.SetRadius(radiusRange[0])
        vtkTube.SetRadiusFactor(radiusRange[0]/float(radiusRange[1]))
        return vtkTube

    def getVtkMapper(self):
        vtkMapper = kwvtk.vtkPolyDataMapper()
        vtkMapper.SetInput(self.getVtkTube().GetOutput())
        vtkMapper.ScalarVisibilityOn()
        vtkMapper.SetScalarModeToUsePointFieldData()
        vtkMapper.ColorByArrayComponent(
            self.getEdgeExtractor().getColorFieldName(),
            self.getEdgeExtractor().getColorFieldComponentIndex(),
        )

        return vtkMapper

    def getActor(self):
        vtkActor = kwvtk.vtkActor()
        vtkActor.SetMapper(self.getVtkMapper())

        return vtkActor

    def addActor(self, addTo):
        addTo.AddActor(self.getActor())


