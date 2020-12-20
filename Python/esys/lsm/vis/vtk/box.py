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

"""Defines the L{Box} and L{Cube} classes"""
from esys.lsm.vis import core
from .modifier     import Modifiable
import vtk as kwvtk

class Box(core.Box,Modifiable):
    def __init__(self, minPt, maxPt):
        core.Box.__init__(self, minPt, maxPt)
        Modifiable.__init__(self)


    def getVtkSource(self):
        box = kwvtk.vtkCubeSource()
        minPt = self.getMinPt()
        maxPt = self.getMaxPt()
        box.SetBounds(
            minPt[0],
            maxPt[0],
            minPt[1],
            maxPt[1],
            minPt[2],
            maxPt[2]
        )

        return box

    def getVtkMapper(self):
        return kwvtk.vtkPolyDataMapper()

    def getActor(self):
        box = self.getVtkSource()

        boxMapper = self.getVtkMapper()
        boxMapper.SetInput(box.GetOutput())

        boxActor = kwvtk.vtkActor()
        boxActor.SetMapper(boxMapper)

        self.applyModifiers(boxActor)

        return boxActor

    def addActor(self, addTo):
        addTo.AddActor(self.getActor())



class Cube(Box):
    """
    Box with all sides the same length.
    """
    def __init__(self, minPt, sideLength):
        """
        Constructs axis-aligned cube.
        @type minPt: sequence of 3 floats
        @param minPt: lower left back corner coordinate.
        @type sideLength: float
        @param sideLength: length of all sides.
        """
        Box.__init__(self, minPt, core.Vec3(minPt)+sideLength)
