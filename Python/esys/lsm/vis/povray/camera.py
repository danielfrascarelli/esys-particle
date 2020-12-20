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

from __future__     import division
from   esys.lsm.vis import core
from   .object       import LightSource

class Camera(core.Camera):
    def __init__(self, posn, lookAt, angle=70.0):
        core.Camera.__init__(self)
        self.posn      = posn
        self.lookAt    = lookAt
        self.angle     = angle
        self.type      = "orthographic"
        self.right     = (-1.0,0.0,0.0) # right-handed coordinate system
        self.up        = ( 0.0,1.0,0.0)
        self.direction = ( 0.0,0.0,0.1)
        self.setZoom(1.0)

    def getPosn(self):
        return self.posn

    def setPosn(self, posn):
        self.posn = posn

    def getLookAt(self):
        return self.lookAt

    def setLookAt(self, lookAt):
        self.lookAt = lookAt

    def setDirection(self, dir):
        self.direction = dir

    def setZoom(self, factor):
        self.zoomFactor = factor

    def getZoom(self):
        return self.zoomFactor

    def setAspect(self, sizeTuple):
        ratio = abs(float(sizeTuple[0])/float(sizeTuple[1]))
        self.right = (-ratio,0.0,0.0) # right handed coordinate system
        self.up    = (0.0,1.0,0.0)

    def getUp(self):
        return core.Vec3(self.up)*(2.05/self.zoomFactor)

    def getRight(self):
        return core.Vec3(self.right)*(2.05/self.zoomFactor)

    def writeSdl(self, f):
        f.write("camera {\n")
        f.write(str(self.type))
        f.write("\nlocation <")
        f.write(str.join(",", (list(map(str,self.posn)))))
        f.write(">\ndirection <")
        f.write(str.join(",", (list(map(str,self.direction)))))
        #f.write(">\nangle ")
        #f.write(str(self.angle))
        f.write(">\nup <")
        f.write(str.join(",", (list(map(str, self.getUp())))))
        f.write(">\nright <")
        f.write(str.join(",", (list(map(str, self.getRight())))))
        f.write(
          ">\nlook_at <"
          +
          str.join(",", (list(map(str,self.lookAt))))
          +
          ">"
        )
        f.write("\n}")

class LightedCamera(Camera):
    def __init__(self, posn, lookAt, angle=70.0):
        Camera.__init__(self, posn, lookAt, angle)
        self.lightSource = LightSource(posn=self.getPosn())

    def setLightSource(self, lightSource):
        self.lightSource = lightSource
        self.lightSource.setPosn(self.getPosn())

    def setPosn(self, posn):
        Camera.setPosn(self, posn)
        self.lightSource.setPosn(posn)

    def writeSdl(self, f):
        Camera.writeSdl(self, f)
        self.lightSource.writeSdl(f)
