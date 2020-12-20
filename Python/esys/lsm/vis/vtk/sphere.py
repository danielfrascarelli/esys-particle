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

from esys.lsm.vis import core
from .modifier     import Modifiable
import vtk as kwvtk

class Sphere(core.Sphere, Modifiable):
    def __init__(self, center, radius):
        core.Sphere.__init__(self, center, radius)
        Modifiable.__init__(self)

    def getVtkSource(self):
        sphere = kwvtk.vtkSphereSource()
        sphere.SetRadius(self.getRadius())
        sphere.SetCenter(self.getCenter())
        sphere.SetPhiResolution(10)
        sphere.SetThetaResolution(10)

        return sphere

    def getVtkMapper(self):
        return kwvtk.vtkPolyDataMapper()

    def getActor(self):
        sphere = self.getVtkSource()

        sphereMapper = self.getVtkMapper()
        sphereMapper.SetInput(sphere.GetOutput())

        sphereActor = kwvtk.vtkActor()
        sphereActor.SetMapper(sphereMapper)

        self.applyModifiers(sphereActor)

        return sphereActor

    def addActor(self, addTo):
        addTo.AddActor(self.getActor())




