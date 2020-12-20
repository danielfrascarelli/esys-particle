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
Defines the L{TriangulatedSurface} class.
"""
from esys.lsm.vis import core
import vtk as kwvtk
from .modifier import Modifiable

class TriangulatedSurface(core.TriangulatedSurface, Modifiable):
    def __init__(self, nodeSequence, faceSequence):
        core.TriangulatedSurface.__init__(self, nodeSequence, faceSequence)
        Modifiable.__init__(self)

    def getVtkIdList(self, idList):
        vtkIdList = kwvtk.vtkIdList()
        for id in idList:
            vtkIdList.InsertNextId(id)
        return vtkIdList

    def getVtkData(self):
        pts = kwvtk.vtkPoints()
        cls = kwvtk.vtkCellArray()

        ptIndex = 0
        for node in self.nodeSequence:
            pts.InsertNextPoint(node)

        for face in self.faceSequence:
            cls.InsertNextCell(self.getVtkIdList(face))

        vtkPolyData = kwvtk.vtkPolyData()
        vtkPolyData.SetPoints(pts)
        vtkPolyData.SetPolys(cls)

        return vtkPolyData

    def getVtkMapper(self):
        mapper = kwvtk.vtkPolyDataMapper()
        mapper.SetInput(self.getVtkData())
        return mapper

    def getActor(self):
        polyMapper = self.getVtkMapper()
        polyActor = kwvtk.vtkActor()
        polyActor.SetMapper(polyMapper)
        self.applyModifiers(polyActor)
        return polyActor

    def addActor(self, addTo):
        addTo.AddActor(self.getActor())
