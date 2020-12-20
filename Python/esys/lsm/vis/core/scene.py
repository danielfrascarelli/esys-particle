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
Defines the L{Scene} base class.
"""
from .exception import raiseNotImplemented
from .color import Colors
from .camera import Camera
from .imageFormat import ImageFormat, PNG, PNM

class RenderDefaults:
    """
    Class for setting default rendering options.
    """
    def __init__(self, defaultOffScreen, defaultSize, defaultInteractive):
        """
        Constructs default values. These values will be returned by the
        L{get} method in for arguments passed with the C{None} value.
        """
        self.defaultOffScreen   = defaultOffScreen
        self.defaultSize        = defaultSize
        self.defaultInteractive = defaultInteractive

    def get(self, offScreen = None, size = None, interactive = None):
        """
        Returns tuple of default options. If an argument is C{None},
        then a default value is returned, if the argument is not C{None}
        then the value is returned unaltered.
        """
        if (offScreen == None):
            offScreen = self.defaultOffScreen
        if (size == None):
            size = self.defaultSize
        if (interactive == None):
            interactive = self.defaultInteractive
        return (offScreen, size, interactive)

_renderDefaults = \
    RenderDefaults(
        defaultOffScreen   = False,
        defaultSize        = [256,256],
        defaultInteractive = True
    )

class Scene(object):
    """
    A container for scene objects.
    """
    def __init__(self, renderDefaults = None):
        """
        Constructs a scene object.
        """
        
        if (renderDefaults == None):
            renderDefaults = _renderDefaults

        self.objectList = []
        self.renderDefaults = renderDefaults

    def initialise(self):
        """
        Initialises the scene, eg sets the background.
        """
        self.setBackground(Colors.White)

    def add(self, object):
        """
        Adds a specified object to the scene.
        @type object: object
        @param object: Object to be rendered in the scene.
        """
        self.objectList.append(object)

    def clear(self):
        """
        Removes all objects from the scene, does not alter camera.
        """
        self.objectList = []

    def setBackground(self, color):
        """
        Sets the background color of the scene.
        @type color: RGB color
        @param color: Backgroung RGB color.
        """
        raiseNotImplemented()

    def getCamera(self):
        """
        Returns a L{Camera} object for this scene.
        @rtype: L{Camera}
        @return: The L{Camera} associated with this scene.
        """

    def getRenderDefaults(
        self,
        offScreen   = None,
        interactive = None,
        size        = None
    ):
        """
        Returns tuple of default render arguments.
        """
        return self.renderDefaults.get(offScreen, interactive, size)

    def render(
        self,
        offScreen   = None,
        interactive = None,
        fileName    = None,
        imageFormat = None,
        size        = None
    ):
        """
        Renders an image of the scene.
        @type offScreen: bool
        @param offScreen: If C{False} the image is displayed on-screen
        in a visible window else the rendered image is not displayed.
        @type interactive: bool
        @param interactive: If C{True} and C{offScreen==True} this method
        will remain I{paused} until a **quit** is received in the image window.
        If C{offScreen==False} this option has no effect.
        @type fileName: str
        @param fileName: The name of the file to which the image is written.
        The C{fileName} extension (eg ".png", ) is used as the image type if
        C{imageFomat} is not specified.
        @type imageFormat: L{ImageFormat}
        @param imageFormat: The type of image produced: L{PNG}, L{PNM}.
        @type size: sequence of 2 int
        @param size: The image width C{size[0]} and height C{size[1]}.
        """
        raiseNotImplemented()
