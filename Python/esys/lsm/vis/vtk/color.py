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

from esys.lsm.vis.core.exception import raiseNotImplemented
from esys.lsm.vis                import core
from .modifier                    import Modifier

class RgbColor(core.RgbColor, Modifier):
    """
    Extends L{core.RgbColor} by implementing the L{applyTo} method.
    """
    def applyTo(self, actor):
        actor.GetProperty().SetColor(self.getRgb())

class Colors:
    def __init__(self):
        raiseNotImplemented(
          "This class is not designed for instantiation."
        )

_nameRgbColorDict = dict()

for attribute in dir(core.Colors):
    attr = getattr(core.Colors, attribute)
    if (isinstance(attr, core.RgbColor)):
        setattr(Colors, attribute, RgbColor(attr[0], attr[1], attr[2]))
        _nameRgbColorDict[str.upper(attr.getName())] = getattr(Colors, attribute)

def findColor(colorName):
    """
    Performs a case-insensitive search for a color with name colorName.
    Returns None if no color matching colorName is found.
    """
    if (str.upper(colorName) in _nameRgbColorDict):
        return _nameRgbColorDict[str.upper(colorName)]
    return None
