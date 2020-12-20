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
#!/bin/env python

from __future__ import division
import sys

from math import sqrt
from optparse import OptionParser
from esys.lsm.util import InstallInfo
if (InstallInfo.haveVtk()):
  from esys.lsm.vis.vtk import *

class Data:
    """
    Objects store position and displacement vectors.
    """
    def __init__(self, point, displacement):
        """
        Initialise data.
        @type point: 3 element sequence of floats
        @param point: position of displacement.
        @type displacement: 3 element sequence of floats
        @param displacement: displacement at point.
        """
        self.point = point
        self.displacement = displacement

    def getPoint(self):
        """
        Returns position.
        @return: data position
        @rtype: 3 element sequence of floats
        """
        return self.point

    def getDisplacement(self):
        """
        Returns displacement.
        @return: data displacement
        @rtype: 3 element sequence of floats
        """
        return self.displacement

    def getDisplacementMagnitude(self):
        """
        Returns displacement magnitude.
        @rtype: float
        @return: magnitude of seld.getDisplacement() vector
        """
        return \
            sqrt(
                sum(
                    [x*x for x in self.displacement]
                )
            )

def surfaceRenderDisplacements(fileName, zField, subsample, scaleMultiplier, imageSize):
    """
    Renders a surface of displacements data.
    @param fileName: name of the data file.
    @type zField: string C{in ["dx", "dy", "dm", "vx", "vy", "vm"]}.
    @param zField: Which component of the data is used for
                  the surface heights.
    @type subsample: int C{>= 1}
    @param subsample: A subset of data is rendered for C{subsample > 1}.
    @type scaleMultiplier: float
    @param scaleMultiplier: Scale factor for surface heights.
    """
    dataList = []
    f = file(fileName, "r")
    xrange=[1000000,-10000000]
    yrange=[1000000,-10000000]
    drange=[1000000,-10000000]
    i = 0
    for line in f.readlines():
        if ((i % subsample) == 0):
            elemList = str.split(str.strip(line), " ")
            elemList = list(map(float, elemList))
            if ((len(zField) > 6) and (zField[0]=="v")):
                dataList.append(
                    Data(elemList[0:3], elemList[6:9])
                )
            else:
                dataList.append(
                    Data(elemList[0:3], elemList[3:6])
                )
        i += 1

    xList = [d.getPoint()[0] for d in dataList]
    yList = [d.getPoint()[1] for d in dataList]
    dList = \
      [d.getDisplacementMagnitude() for d in dataList] \
      +\
      [d.getDisplacement()[0] for d in dataList] \
      +\
      [d.getDisplacement()[1] for d in dataList]

    xrange = [min(xList), max(xList)]
    yrange = [min(yList), max(yList)]
    drange = [min(dList), max(dList)]

    scaleFactor = 1.0
    if (drange[1] > 0.0):
        scaleFactor =       \
            scaleMultiplier \
            *               \
            (max([xrange[1]-xrange[0], yrange[1]-yrange[0]])/2.0)/drange[1]

    if (zField[1] == "x"):
        extr = lambda x : [x.getPoint()[0], x.getPoint()[1], scaleFactor*x.getDisplacement()[0]]
    elif (zField[1] == "y"):
        extr = lambda x : [x.getPoint()[0], x.getPoint()[1], scaleFactor*x.getDisplacement()[1]]
    else:
        extr = lambda x : [x.getPoint()[0], x.getPoint()[1], scaleFactor*x.getDisplacementMagnitude()]
    
    scene = Scene()
    scene.setBackground(Colors.LightBlue)
    scene.add(SurfaceData(dataList, PointExtractor(extr)))
    scene.render(offScreen=False, size=imageSize)

def getOptionParser():
    """
    Returns an object for parsing command line arguemnts.
    @rtype: C{optparse.OptionParser}
    @return: a parser initialised with some relevent command line options.
    """
    usage = "usage: %prog [options] <dataFile>\n\n"
    usage += "Displays a fitted surface of displacement/velocity data"
    usage += " from the\nfile"
    usage += " <dataFile>. This file is expected to contain lines\nof"
    usage += " the form \"px py pz dx dy dz vx vy vz\"."
    usage += " The pz and vz components\nare ignored."
    parser = OptionParser(usage=usage)
    parser.add_option(
      "-z", "--z-component",
      dest="zComponent",
      metavar="Z",
      default="dmag",
      help=\
          "Defines the component of the displacement or velocity which " +\
          "is used as the surface height," +\
          "Z={'dx', 'dy', 'dmagnitude', 'vx', 'vy', 'vmagnitude'}" +\
          " (default Z=\"dmagnitude\")."
    )
    parser.add_option(
      "-s", "--subsample",
      dest="subsample",
      metavar="S",
      type = "int",
      default=1,
      help=\
          "Only every S-th data point is read from file, rest are skipped " +\
          " (default S=1)."
    )
    parser.add_option(
      "-x", "--scale-multiplier",
      dest="scaleMultiplier",
      metavar="S",
      type = "float",
      default=1.0,
      help=\
          "Scale heights by this amount." +\
          " (default S=1.0)."
    )
    parser.add_option(
      "-i", "--image-size",
      dest="imageSize",
      metavar="P",
      type = "int",
      default=1024,
      help=\
          "Size of the image will be PxP pixels." +\
          " (default P=1024)."
    )

    return parser

__doc__ = \
"""
Defines functions for rendering displacement data as a 3D surface.
This module can be as a C{{__main__}} and has the following usage::
{0:s}
""".format("  " + str.replace(getOptionParser().format_help(), "\n", "\n  "))

if (__name__ == "__main__"):
    parser = getOptionParser()
    (options, args) = parser.parse_args(sys.argv[1:])
   
    fileName = args[0]
    surfaceRenderDisplacements(
        fileName,
        options.zComponent,
        options.subsample,
        options.scaleMultiplier,
        [options.imageSize, options.imageSize]
    )
