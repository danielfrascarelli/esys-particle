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

from __future__ import print_function

import string
import os
import os.path
import itertools

from   esys.lsm.vis import core
from   .povRenderer  import PovRenderer
from   .object       import Background
from   .camera       import LightedCamera
from   .color        import Colors, RgbColor

class Renderer:
    def __init__(self):
        self.povRenderer = None
        self.background = Background(color=Colors.White)
        self.camera = self.getDefaultCamera()

    def getDefaultCamera(self):
        return LightedCamera(posn=[0,0,-10], lookAt=[0,0,0])

    def getCamera(self):
        return self.camera

    def setCamera(self, camera):
        self.camera = camera

    def getPovRenderer(self):
        if (self.povRenderer == None):
            self.povRenderer = PovRenderer()
        return self.povRenderer

    def setOffScreen(self, doOffScreen):
        self.getPovRenderer().setOffScreen(doOffScreen)

    def getOffScreen(self):
        return self.getPovRenderer().getOffScreen()

    def setInteractive(self, interactive):
        self.getPovRenderer().setInteractive(interactive)

    def getInteractive(self):
        return self.getPovRenderer().getInteractive()

    def setSize(self, size):
        self.getPovRenderer().setSize(size)

    def getSize(self):
        return self.getPovRenderer().getSize()

    def clear(self):
        pass

    def setBackground(self, color):
        self.background = Background(color)

    def getBackground(self):
        return self.background

    def getRendererObjectList(self):
        return [self.getBackground(), self.getCamera()]

    def getImageFormat(self, fileName):
        return core.getFormatFromExtension(fileName)

    def getObjectList(self):
        return self.objectList

    def setObjectList(self, objectList):
        self.objectList = objectList

    def writeSdl(self, f, objectList):
        self.getCamera().setAspect(self.getSize())
        for o in itertools.chain(self.getRendererObjectList(), objectList):
          o.writeSdl(f)
          print("\n",file=f)
      
    def render(self, objectList, fileName = None, imageFormat = None):
        if (fileName != None):
            if (imageFormat == None):
                imageFormat = self.getImageFormat(fileName)
        self.getPovRenderer().setFileName(fileName, imageFormat)
        f = self.getPovRenderer().openInput()
        self.writeSdl(f, objectList)
        self.getPovRenderer().closeInput()

class RenderDefaults(core.RenderDefaults):
    pass

class Scene(core.Scene):
    def __init__(self, renderDefaults=None):
        self.renderer = None
        core.Scene.__init__(self, renderDefaults)
        self.initialise()

    def getRenderer(self):
        if (self.renderer == None):
            self.renderer = Renderer()
        return self.renderer

    def setBackground(self, color):
        self.getRenderer().setBackground(
          RgbColor(color[0],color[1],color[2],color.getName())
        )

    def getCamera(self):
        return self.getRenderer().getCamera()

    def setCamera(self, camera):
        return self.getRenderer().setCamera(camera)

    def clear(self):
        core.Scene.clear(self)
        self.getRenderer().clear()

    def render(
        self,
        offScreen   = None,
        interactive = None,
        fileName    = None,
        imageFormat = None,
        size        = None
    ):
        (offScreen, size, interactive) = \
            self.getRenderDefaults(offScreen, size, interactive)
        self.getRenderer().clear()
        self.getRenderer().setOffScreen(offScreen)
        self.getRenderer().setInteractive(interactive)
        self.getRenderer().setSize(size)
        self.getRenderer().setObjectList(self.objectList)
        self.getRenderer().render(self.objectList, fileName, imageFormat)
