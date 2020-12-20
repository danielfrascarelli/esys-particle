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
Defines L{RgbColor} base class and L{Colors} pseudo class.
"""
from .exception import raiseNotImplemented

class RgbColor(object):
    """
    Represents a named color as a C{(red, green, blue)} tuple.
    """
    def __init__(self, red, green, blue, name=None):
        """
        Initialises a color object.
        @type red: float
        @param red: Red value between 0.0 and 1.0 inclusive.
        @type green: float
        @param green: Green value between 0.0 and 1.0 inclusive.
        @type blue: float
        @param blue: Blue value between 0.0 and 1.0 inclusive.
        @type name: str
        @param name: The name assigned to this color object.
        """
        self.rgb  = (red, green, blue)
        self.name = name

    def __len__(self):
        return 3

    def __getitem__(self, i):
        return self.rgb[i]

    def getRgb(self):
        """
        Returns red, green and blue intensities as a C{(r,g,b)} tuple.
        @rtype: tuple
        @return: C{(r,g,b)} tuple.
        """
        return self.rgb

    def getName(self):
        """
        Returns the name associated with this color.
        @rtype: str
        @return: The name associated with this color.
        """
        return self.name

    def setName(self, name):
        """
        Sets the name associated with this color.
        @type name: str
        @param name: Name assigned to this color.
        """
        self.name = name

    def __mul__(self, fraction):
        return \
            RgbColor(
                self.rgb[0]*fraction,
                self.rgb[1]*fraction,
                self.rgb[2]*fraction
            )

class Colors:
    """
    Defines class-scope variables representing particular
    color instances.
    """
    def __init__(self):
        """
        Raises C{NotImplementedError}, this class may not
        be instantiated.
        @raise NotImplementedError: Always.
        """
        raiseNotImplementedError(
          "This class is not designed for instantiation."
        )
    #
    # From Povray colors.inc
    #
    Red     = RgbColor(1, 0, 0)
    Green   = RgbColor(0, 1, 0)
    Blue    = RgbColor(0, 0, 1)
    Yellow  = RgbColor(1, 1, 0)
    Cyan    = RgbColor(0, 1, 1)
    Magenta = RgbColor(1, 0, 1)
    
    White   = RgbColor(1, 1, 1)
    Black   = RgbColor(0, 0, 0)
    
    # These grays are useful for fine-tuning lighting RgbColor(values
    # and for other areas where subtle variations of grays are needed.
    # PERCENTAGE GRAYS:
    Gray05 = White*0.05
    Gray10 = White*0.10
    Gray15 = White*0.15
    Gray20 = White*0.20
    Gray25 = White*0.25
    Gray30 = White*0.30
    Gray35 = White*0.35
    Gray40 = White*0.40
    Gray45 = White*0.45
    Gray50 = White*0.50
    Gray55 = White*0.55
    Gray60 = White*0.60
    Gray65 = White*0.65
    Gray70 = White*0.70
    Gray75 = White*0.75
    Gray80 = White*0.80
    Gray85 = White*0.85
    Gray90 = White*0.90
    Gray95 = White*0.95
    
    # OTHER GRAYS
    DimGray = RgbColor(0.329412, 0.329412, 0.329412)
    DimGrey = RgbColor(0.329412, 0.329412, 0.329412)
    Gray = RgbColor(0.752941, 0.752941, 0.752941)
    Grey = RgbColor(0.752941, 0.752941, 0.752941)
    LightGray = RgbColor(0.658824, 0.658824, 0.658824)
    LightGrey = RgbColor(0.658824, 0.658824, 0.658824)
    VLightGray = RgbColor(0.80, 0.80, 0.80)
    VLightGrey = RgbColor(0.80, 0.80, 0.80)
    
    Aquamarine = RgbColor(0.439216, 0.858824, 0.576471)
    BlueViolet = RgbColor(0.62352, 0.372549, 0.623529)
    Brown = RgbColor(0.647059, 0.164706, 0.164706)
    CadetBlue = RgbColor(0.372549, 0.623529, 0.623529)
    Coral = RgbColor(1.0, 0.498039, 0.0)
    CornflowerBlue = RgbColor(0.258824, 0.258824, 0.435294)
    DarkGreen = RgbColor(0.184314, 0.309804, 0.184314)
    DarkOliveGreen = RgbColor(0.309804, 0.309804, 0.184314)
    DarkOrchid = RgbColor(0.6, 0.196078, 0.8)
    DarkSlateBlue = RgbColor(0.119608, 0.137255, 0.556863)
    DarkSlateGray = RgbColor(0.184314, 0.309804, 0.309804)
    DarkSlateGrey = RgbColor(0.184314, 0.309804, 0.309804)
    DarkTurquoise = RgbColor(0.439216, 0.576471, 0.858824)
    Firebrick = RgbColor(0.556863, 0.137255, 0.137255)
    ForestGreen = RgbColor(0.137255, 0.556863, 0.137255)
    Gold = RgbColor(0.8, 0.498039, 0.196078)
    Goldenrod = RgbColor(0.858824, 0.858824, 0.439216)
    GreenYellow = RgbColor(0.576471, 0.858824, 0.439216)
    IndianRed = RgbColor(0.309804, 0.184314, 0.184314)
    Khaki = RgbColor(0.623529, 0.623529, 0.372549)
    LightBlue = RgbColor(0.74902, 0.847059, 0.847059)
    LightSteelBlue = RgbColor(0.560784, 0.560784, 0.737255)
    LimeGreen = RgbColor(0.196078, 0.8, 0.196078)
    Maroon = RgbColor(0.556863, 0.137255, 0.419608)
    MediumAquamarine = RgbColor(0.196078, 0.8, 0.6)
    MediumBlue = RgbColor(0.196078, 0.196078, 0.8)
    MediumForestGreen = RgbColor(0.419608, 0.556863, 0.137255)
    MediumGoldenrod = RgbColor(0.917647, 0.917647, 0.678431)
    MediumOrchid = RgbColor(0.576471, 0.439216, 0.858824)
    MediumSeaGreen = RgbColor(0.258824, 0.435294, 0.258824)
    MediumSlateBlue = RgbColor(0.498039, 0.0, 1.0)
    MediumSpringGreen = RgbColor(0.498039, 1.0, 0.0)
    MediumTurquoise = RgbColor(0.439216, 0.858824, 0.858824)
    MediumVioletRed = RgbColor(0.858824, 0.439216, 0.576471)
    MidnightBlue = RgbColor(0.184314, 0.184314, 0.309804)
    Navy = RgbColor(0.137255, 0.137255, 0.556863)
    NavyBlue = RgbColor(0.137255, 0.137255, 0.556863)
    Orange = RgbColor(1, 0.5, 0.0)
    OrangeRed = RgbColor(1.0, 0.25, 0.0)
    Orchid = RgbColor(0.858824, 0.439216, 0.858824)
    PaleGreen = RgbColor(0.560784, 0.737255, 0.560784)
    Pink = RgbColor(0.737255, 0.560784, 0.560784)
    Plum = RgbColor(0.917647, 0.678431, 0.917647)
    Salmon = RgbColor(0.435294, 0.258824, 0.258824)
    SeaGreen = RgbColor(0.137255, 0.556863, 0.419608)
    Sienna = RgbColor(0.556863, 0.419608, 0.137255)
    SkyBlue = RgbColor(0.196078, 0.6, 0.8)
    SlateBlue = RgbColor(0.0, 0.498039, 1.0)
    SpringGreen = RgbColor(0.0, 1.0, 0.498039)
    SteelBlue = RgbColor(0.137255, 0.419608, 0.556863)
    Tan = RgbColor(0.858824, 0.576471, 0.439216)
    Thistle = RgbColor(0.847059, 0.74902, 0.847059)
    Turquoise = RgbColor(0.678431, 0.917647, 0.917647)
    Violet = RgbColor(0.309804, 0.184314, 0.309804)
    VioletRed = RgbColor(0.8, 0.196078, 0.6)
    Wheat = RgbColor(0.847059, 0.847059, 0.74902)
    YellowGreen = RgbColor(0.6, 0.8, 0.196078)
    SummerSky = RgbColor(0.22, 0.69, 0.87)
    RichBlue = RgbColor(0.35, 0.35, 0.67)
    Brass =  RgbColor(0.71, 0.65, 0.26)
    Copper = RgbColor(0.72, 0.45, 0.20)
    Bronze = RgbColor(0.55, 0.47, 0.14)
    Bronze2 = RgbColor(0.65, 0.49, 0.24)
    Silver = RgbColor(0.90, 0.91, 0.98)
    BrightGold = RgbColor(0.85, 0.85, 0.10)
    OldGold =  RgbColor(0.81, 0.71, 0.23)
    Feldspar = RgbColor(0.82, 0.57, 0.46)
    Quartz = RgbColor(0.85, 0.85, 0.95)
    Mica = Black
    NeonPink = RgbColor(1.00, 0.43, 0.78)
    DarkPurple = RgbColor(0.53, 0.12, 0.47)
    NeonBlue = RgbColor(0.30, 0.30, 1.00)
    CoolCopper = RgbColor(0.85, 0.53, 0.10)
    MandarinOrange = RgbColor(0.89, 0.47, 0.20)
    LightWood = RgbColor(0.91, 0.76, 0.65)
    MediumWood = RgbColor(0.65, 0.50, 0.39)
    DarkWood = RgbColor(0.52, 0.37, 0.26)
    SpicyPink = RgbColor(1.00, 0.11, 0.68)
    SemiSweetChoc = RgbColor(0.42, 0.26, 0.15)
    BakersChoc = RgbColor(0.36, 0.20, 0.09)
    Flesh = RgbColor(0.96, 0.80, 0.69)
    NewTan = RgbColor(0.92, 0.78, 0.62)
    NewMidnightBlue = RgbColor(0.00, 0.00, 0.61)
    VeryDarkBrown = RgbColor(0.35, 0.16, 0.14)
    DarkBrown = RgbColor(0.36, 0.25, 0.20)
    DarkTan = RgbColor(0.59, 0.41, 0.31)
    GreenCopper = RgbColor(0.32, 0.49, 0.46)
    DkGreenCopper = RgbColor(0.29, 0.46, 0.43)
    DustyRose = RgbColor(0.52, 0.39, 0.39)
    HuntersGreen = RgbColor(0.13, 0.37, 0.31)
    Scarlet = RgbColor(0.55, 0.09, 0.09)
    
    MedPurple =  RgbColor(0.73, 0.16, 0.96)
    LightPurple = RgbColor(0.87, 0.58, 0.98)
    VeryLightPurple = RgbColor(0.94, 0.81, 0.99)

_nameRgbColorDict = dict()

for attribute in dir(Colors):
    attr = getattr(Colors, attribute)
    if (isinstance(attr, RgbColor)):
        attr.setName(attribute)
        _nameRgbColorDict[str.upper(attr.getName())] = attr

def findColor(colorName):
    """
    Performs a case-insensitive search for a color with name colorName.
    Returns None if no color matching colorName is found.
    """
    if (str.upper(colorName) in _nameRgbColorDict):
        return _nameRgbColorDict[str.upper(colorName)]
    return None
