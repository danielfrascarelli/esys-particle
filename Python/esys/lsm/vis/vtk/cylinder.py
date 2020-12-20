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
from .util         import computeRotateWxyz
from .modifier     import Modifiable
import vtk as kwvtk
import math

class Cylinder(core.Cylinder, Modifiable):
    def __init__(self, endPt1, endPt2, radius):
        core.Cylinder.__init__(self, endPt1, endPt2, radius)
        Modifiable.__init__(self)

    def getVtkSource(self):
        cylinder = kwvtk.vtkCylinderSource()
        cylinder.SetRadius(self.getRadius())
        cylinder.SetCenter(self.getCenter())
        cylinder.SetHeight(self.getHeight())
        cylinder.SetResolution(10)

        return cylinder

    def getVtkMapper(self):
        return kwvtk.vtkPolyDataMapper()

    def getActor(self):
        cylinder = self.getVtkSource()

        cylinderMapper = self.getVtkMapper()
        cylinderMapper.SetInput(cylinder.GetOutput())

        cylinderActor = kwvtk.vtkActor()
        cylinderActor.SetMapper( cylinderMapper )
        cylinderActor.SetOrigin(self.getCenter())

        yDir = core.Vec3(0,1,0) # original vtk cylinder-source axis direction
        axisDir = core.Vec3(self.getEndPt2())-core.Vec3(self.getEndPt1())
        (theta, rotAxis) = computeRotateWxyz(yDir, axisDir)
        cylinderActor.RotateWXYZ(
            theta,
            rotAxis[0],
            rotAxis[1],
            rotAxis[2]
        )
        self.applyModifiers(cylinderActor)

        return cylinderActor

    def addActor(self, addTo):
        addTo.AddActor(self.getActor())

class Disk(Cylinder):
    """
    Objects of this class represent cylinders.
    """
    def __init__(self, center, radius, height, direction=None):
        if (direction == None):
            direction = core.Vec3(0,0,1)
        direction = core.Vec3(direction)
        nrm = direction.norm()
        if (nrm == 0.0):
            direction = core.Vec3(0,0,1)
            nrm = 1.0

        center = core.Vec3(center)
        direction = direction/nrm
        Cylinder.__init__(
            self,
            endPt1 = center - (direction*(height*0.5)),
            endPt2 = center + (direction*(height*0.5)),
            radius = radius
        )
