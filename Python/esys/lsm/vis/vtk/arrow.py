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

"""Defines the L{Arrow} class"""
from __future__   import division
from esys.lsm.vis import core
from .util         import computeRotateWxyz
from .modifier     import Modifiable
import vtk as kwvtk # Kitware's vtk package

class Arrow(core.Arrow, Modifiable):
    def __init__(self, tailPt, headPt, tailRadius, headRadius, headLength):
        core.Arrow.__init__(
            self,
            tailPt,
            headPt,
            tailRadius,
            headRadius,
            headLength
        )
        Modifiable.__init__(self)

    def getVtkSource(self):
        arrow = kwvtk.vtkArrowSource()
        scale = self.getLength()
        arrow.SetTipLength(self.getHeadLength()/scale)
        arrow.SetTipRadius(self.getHeadRadius()/scale)
        arrow.SetShaftRadius(self.getTailRadius()/scale)
        arrow.SetTipResolution(10)
        arrow.SetShaftResolution(10)

        return arrow

    def getVtkMapper(self):
        return kwvtk.vtkPolyDataMapper()

    def getActor(self):
        arrow = self.getVtkSource()

        arrowMapper = self.getVtkMapper()
        arrowMapper.SetInput(arrow.GetOutput())

        arrowActor = kwvtk.vtkActor()
        arrowActor.SetMapper(arrowMapper)
        arrowLength = self.getLength()
        arrowActor.SetScale([arrowLength]*3)

        (theta, axis) = \
            computeRotateWxyz(
                (1,0,0),
                core.Vec3(self.getHeadPt())-core.Vec3(self.getTailPt())
            )
        arrowActor.RotateWXYZ(theta, axis[0], axis[1], axis[2])
        arrowActor.SetPosition(self.getTailPt())

        self.applyModifiers(arrowActor)

        return arrowActor

    def addActor(self, addTo):
        addTo.AddActor(self.getActor())









