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

from .modifier import Pigment
from .color    import Colors

class EdgeExtractor(core.EdgeExtractor):
    def __init__(
        self,
        pointListMap    = lambda dataRecord: dataRecord.getPointList(),
        radiusListMap   = lambda dataRecord: 1.0,
        colorValListMap = lambda dataRecord: Pigment(Colors.Red),
        radiusScale     = 1.0
    ):
        core.EdgeExtractor.__init__(
            self,
            pointListMap,
            radiusListMap,
            colorValListMap,
            radiusScale
        )

    def getResizedList(self, lst, newLen):
        if (not hasattr(lst, "__iter__")):
            newLst = [lst]
        else:
            newLst = list(lst)
        if (len(newLst) < newLen):
            newLst = newLst*((newLen//(len(newLst))) + 1)
        del newLst[newLen:]
        return newLst

    def writeSdl(self, f, record):
        ptList = self.getPointList(record)
        numPts = len(ptList)
        radiusList = \
          self.getResizedList(self.getRadiusList(record), numPts)
        colorList  = \
          self.getResizedList(self.getColorValList(record), numPts)
        for i in range(0, numPts-1):
            f.write("cone {\n<")
            f.write(str.join(",", (list(map(str,ptList[i])))))
            f.write(">,")
            f.write(str(radiusList[i]*self.getRadiusScale()))
            f.write("\n<")
            f.write(str.join(",", (list(map(str,ptList[i+1])))))
            f.write(">,")
            f.write(str(radiusList[i+1]*self.getRadiusScale()))
            f.write("\n")
            if (hasattr(colorList[i], "writeSdl")):
                colorList[i].writeSdl(f)
            f.write("}\n")
