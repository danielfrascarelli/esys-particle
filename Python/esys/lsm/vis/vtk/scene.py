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
"""
import string
import os
import os.path

from esys.lsm.vis import core
from .camera       import Camera

import vtk as kwvtk

formatWriterDict = dict()
formatWriterDict[str(core.PNG)] = kwvtk.vtkPNGWriter
formatWriterDict[str(core.PNM)] = kwvtk.vtkPNMWriter

class Renderer:
    def __init__(self):
        self.vtkRenderWindow = None
        self.vtkRenderer = None
        self.vtkInteractor = None
        self.camera = None
        self.interactive = True

    def getVtkRenderWindow(self):
        if (self.vtkRenderWindow == None):
            self.vtkRenderWindow = self.createVtkRenderWindow()
            self.vtkRenderWindow.AddRenderer(self.getVtkRenderer())
        return self.vtkRenderWindow

    def createVtkRenderWindow(self):
        return kwvtk.vtkRenderWindow()

    def createVtkRenderer(self):
        return kwvtk.vtkRenderer()

    def getVtkRenderer(self):
        if (self.vtkRenderer == None):
            self.vtkRenderer = self.createVtkRenderer()
        return self.vtkRenderer

    def setOffScreen(self, doOffScreen):
        self.getVtkRenderWindow().SetOffScreenRendering(int(doOffScreen))

    def getOffScreen(self):
        return \
            bool(
              self.getVtkRenderWindow().GetOffScreenRendering()
            )

    def getInteractive(self):
        return self.interactive

    def setInteractive(self, interactive):
        self.interactive = bool(interactive)

    def getCamera(self):
        if (self.camera == None):
            self.camera = Camera(self.getVtkRenderer().GetActiveCamera())
        return self.camera

    def getVtkInteractor(self):
        if (self.vtkInteractor == None):
            self.vtkInteractor = kwvtk.vtkRenderWindowInteractor()
            self.vtkInteractor.SetInteractorStyle(
                kwvtk.vtkInteractorStyleTrackballCamera()
            )
            

        return self.vtkInteractor

    def updateInteractor(self):
        if ((not self.getOffScreen()) and self.getInteractive()):
            if (not self.getVtkInteractor().GetRenderWindow()):
              self.getVtkInteractor().SetRenderWindow(self.getVtkRenderWindow())
            if (not self.getVtkInteractor().GetInitialized()):
                self.getVtkInteractor().Initialize()
            if (not self.getVtkInteractor().GetEnabled()):
                self.getVtkInteractor().Enable()
        elif (self.getVtkInteractor().GetEnabled()):
            self.getVtkInteractor().Disable()

    def getOffScreen(self):
        return self.getVtkRenderWindow().GetOffScreenRendering()

    def setSize(self, size):
        if (size != None):
            self.getVtkRenderWindow().SetSize(size)

    def getSize(self):
        return self.getVtkRenderWindow().GetSize()

    def getImageWriter(self, imageFormat):
        return formatWriterDict[str.upper(str(imageFormat))]()

    def getImageFormat(self, fileName):
        return core.getFormatFromExtension(fileName)

    def clear(self):
        self.getVtkRenderer().RemoveAllProps()
        self.getVtkRenderer().GetProps().RemoveAllItems()

    def add(self, iteratable):
        vtkRenderer = self.getVtkRenderer()
        for thing in iteratable:
            thing.addActor(vtkRenderer)

    def setBackground(self, color):
        self.getVtkRenderer().SetBackground(color.getRgb())

    def startInteractor(self):
        if (self.getVtkInteractor().GetEnabled()):
            self.getVtkInteractor().Start()

    def render(self, fileName = None, imageFormat = None):
        self.getVtkRenderWindow().Render()
        if (fileName != None):
            if (imageFormat == None):
                imageFormat = self.getImageFormat(fileName)
            imageInput = kwvtk.vtkRenderLargeImage()
            imageInput.SetMagnification(1)
            imageInput.SetInput(self.getVtkRenderer())
            writer = self.getImageWriter(imageFormat)
            writer.SetInput(imageInput.GetOutput())
            writer.SetFileName(fileName)
            writer.Write()

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
        self.getRenderer().setBackground(color)

    def clear(self):
        core.Scene.clear(self)
        self.getRenderer().clear()

    def getCamera(self):
        return self.getRenderer().getCamera()

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
        self.getRenderer().add(self.objectList)
        self.getRenderer().setOffScreen(offScreen)
        self.getRenderer().setInteractive(interactive)
        self.getRenderer().setSize(size)
        self.getRenderer().updateInteractor()
        self.getRenderer().render(fileName, imageFormat)
        self.getRenderer().startInteractor()
