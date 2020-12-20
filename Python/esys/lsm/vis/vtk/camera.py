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

class Camera(core.Camera):
    def __init__(self, vtkCamera):
        core.Camera.__init__(self)
        self.vtkCamera = vtkCamera
        self.SetUseHorizontalViewAngle(True)
        self.SetViewAngle(70.0)
        self.SetParallelProjection(True)
        self.SetClippingRange(0.0001, 1.0e8)
        self.setZoom(1.0)

    def __getattr__(self, attr):
        return getattr(self.vtkCamera, attr)

    def getPosn(self):
        return self.vtkCamera.GetPosition()

    def setZoom(self, factor):
        self.zoomFactor = factor
        self.Zoom(factor)

    def getZoom(self):
        return self.zoomFactor

    def setPosn(self, posn):
        self.vtkCamera.SetPosition(posn)

    def getLookAt(self):
        return self.vtkCamera.GetFocalPoint()

    def setLookAt(self, lookAt):
        self.vtkCamera.SetFocalPoint(lookAt)

    def setAngle(self, angle):
        self.vtkCamera.SetViewAngle(angle)

    def getAngle(self):
        return self.vtkCamera.GetViewAngle()
