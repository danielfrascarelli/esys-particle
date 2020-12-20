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
Defines the L{ImageFormat} class and functions for mapping a
file name extension to an associated C{ImageFormat} object.
"""
import os
import os.path

class ImageFormat(object):
    """
    Class representing an image format.
    """
    def __init__(self, name):
        """
        Constructor.
        @type name: str
        @param name: Name assigned to this image format.
        """
        self.name = name

    def getName(self):
        """
        Returns the name associated with this image format.
        @rtype: str
        """
        return self.name

    def __str__(self):
        return self.getName()


PNG = ImageFormat("PNG")
PNM = ImageFormat("PNM")

_nameFormatDict = dict()
_nameFormatDict[str.upper(str(PNG))] = PNG
_nameFormatDict[str.upper(str(PNM))] = PNM

def _getDelimitedFormatNameString():
    return ", ".join(map(str,list(_nameFormatDict.keys())))

def getFormatFromName(formatName, ext=None):
    """
    Returns the C{{ImageFormat}} object which corresponds
    to a specified image-format name (string).
    @type formatName: str
    @param formatName: The name of an image format, one of: {0:s}
    @type ext: str
    @param ext: File name extension for error message string.
    """.format(_getDelimitedFormatNameString())
    if str.upper(formatName) in _nameFormatDict:
        return _nameFormatDict[str.upper(formatName)]
    raise \
      ValueError(
        (
            "No image format found which matched extension '{0:s}';" +
            " valid image file formats are: {1:s}"
        ).format(ext, _getDelimitedFormatNameString())
      )

def getFormatFromExtension(fileName):
    """
    Returns the C{ImageFormat} object which corresponds
    to a specified file name. Uses the C{fileName} extension
    to try and deduce the corresponding C{ImageFormat} object.
    @type fileName: str
    @param fileName: A file name.
    @rtype: C{ImageFormat}
    @return: An C{ImageFormat} object corresponding to the
    specified file name (and corresponding file name extension).
    """
    (base, ext) = os.path.splitext(fileName)
    if (len(ext) > 0):
        formatName = str.lstrip(ext, ".")
    else:
        raise ValueError(
            "Could not determine image format from file "
            +
            "name " + fileName + ", no extension."
        )
    return getFormatFromName(formatName, ext)
