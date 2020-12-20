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
Defines classes, instances of which, can be used to generate
spherical grains packed with randomly sized circular particles.
"""
from esys.lsm                     import Logging
from esys.lsm.util                import *
from esys.lsm.geometry.GeometryPy import *

def getLogger():
    return Logging.getLogger("esys.lsm.geometry")

class ScaledMinMaxGrainGen(RndGrainGen):
    """
    Grain generator which can be used in conjunction with
    the L{GougeConfigPrms} and L{GrainRandomBoxPacker} classes.
    Generates grains which consist of a particles packed into
    a circle/sphere. The minimum and maximum particle radii
    are scaled with respect to grain size.
    @see: getGrain
    """
    def __init__(
        self,
        minGRadius,
        maxGRadius,
        do2d=False,
        radiusScale=(0.1, 0.4),
        maxInsertFails=5000
    ):
        """
        Initialise the generator.
        @type minGRadius: float
        @param minGRadius: minimum grain radius
        @type maxGRadius: float
        @param maxGRadius: maximum grain radius
        @type do2d: bool
        @param do2d: If True, grain particles all occur
          on the M{z=0} plane (ie particles packed into
          a circle).
        @type radiusScale: tuple of two floats
        @param radiusScale: Randomly generated spherical
          particles will have radii which will lie in the range
          C{(radiusScale[0]*minGRadius,radiusScale[1]*maxGRadius)}.
        """
        self.do2d = do2d
        self.maxInsertFails = maxInsertFails
        self.minPRadiusScale = radiusScale[0]
        self.maxPRadiusScale = radiusScale[1]
        minPRadius = self.minPRadiusScale*minGRadius
        maxPRadius = self.maxPRadiusScale*maxGRadius
        RndGrainGen.__init__(
            self,
            minGRadius,
            maxGRadius,
            minPRadius,
            maxPRadius
        )

    def getGrain(self, p):
        """
        Uses L{RandomSpherePacker} to create a grain consisting
        of randomly sized particles packed into a sphere with
        centre C{p.getPosn()} and radius C{p.getRadius()}. The
        random particle radii will lie in the range
        C{(self.minPRadiusScale*p.getRadius(),self.maxPRadiusScale*p.getRadius())}
        @type p: L{SimpleSphere}
        @param p: C{p.getPosn()} and C{p.getRadius()} determine
          the size and position of the returned grain.
        @rtype: L{Grain}
        @return: A I{spherical} grain of randomly sized particles.
        """
        minRadius=p.getRadius()*self.minPRadiusScale
        packer = \
          RandomSpherePacker(
            minRadius=minRadius,
            maxRadius=p.getRadius()*self.maxPRadiusScale,
            cubicPackRadius=p.getRadius()*self.maxPRadiusScale,
            maxInsertFails=self.maxInsertFails,
            bSphere=BoundingSphere(p.getPosn(), p.getRadius()),
            do2d=self.do2d,
            tolerance = 0.1*minRadius
          )
        packer.generate()
        g = Grain()
        for p in packer.getParticleIterator():
            g.createParticle(p)
        getLogger().info(
          "Created grain of " + str(g.getNumParticles()) + " particles"
        )
        return g

class ScaledMaxGrainGen(RndGrainGen):
    """
    Grain generator which can be used in conjunction with
    the L{GougeConfigPrms} and L{GrainRandomBoxPacker} classes.
    Generates grains which consist of a particles packed into
    a circle/sphere. The maximum particle radius
    is scaled with respect to grain size, the minimum
    particle radius is independent of grain size.
    @see: getGrain
    """
    def __init__(
        self,
        minGRadius,
        maxGRadius,
        do2d=False,
        minPRadius=None,
        radiusScale=0.1,
        maxInsertFails=5000
    ):
        """
        Initialise the generator.
        @type minGRadius: float
        @param minGRadius: minimum grain radius
        @type maxGRadius: float
        @param maxGRadius: maximum grain radius
        @type do2d: bool
        @param do2d: If True, grain particles all occur
          on the M{z=0} plane (ie particles packed into
          a circle).
        @type minPRadius: float
        @param minPRadius: The minimum radius of randomly
          generated particles.
        @type radiusScale: float
        @param radiusScale: Randomly generated spherical
          particles will have radii which will lie in the range
          C{(minPRadius,radiusScale*maxGRadius)}.
        """
        if (minPRadius == None):
            minPRadius = minGRadius*0.15
        self.do2d            = do2d
        self.maxInsertFails  = maxInsertFails
        self.minPRadius      = minPRadius
        self.maxPRadiusScale = radiusScale
        maxPRadius = self.maxPRadiusScale*maxGRadius
        RndGrainGen.__init__(
            self,
            minGRadius,
            maxGRadius,
            minPRadius,
            maxPRadius
        )

    def getGrain(self, p):
        """
        Uses L{RandomSpherePacker} to create a grain consisting
        of randomly sized particles packed into a sphere with
        centre C{p.getPosn()} and radius C{p.getRadius()}. The
        random particle radii will lie in the range
        C{(self.minPRadius,self.maxPRadiusScale*p.getRadius())}
        @type p: L{SimpleSphere}
        @param p: C{p.getPosn()} and C{p.getRadius()} determine
          the size and position of the returned grain.
        @rtype: L{Grain}
        @return: A I{spherical} grain of randomly sized particles.
        """
        minRadius = self.minPRadius
        maxRadius = p.getRadius()*self.maxPRadiusScale
        packer = \
          RandomSpherePacker(
            minRadius=minRadius,
            maxRadius=maxRadius,
            cubicPackRadius=maxRadius,
            maxInsertFails=self.maxInsertFails,
            bSphere=BoundingSphere(p.getPosn(), p.getRadius()),
            do2d=self.do2d,
            tolerance = 0.1*minRadius
          )
        packer.generate()
        g = Grain()
        for p in packer.getParticleIterator():
            g.createParticle(p)
        getLogger().info(
          "Created grain of " + str(g.getNumParticles()) + " particles"
        )
        return g

class AbsGrainGen(RndGrainGen):
    """
    Grain generator which can be used in conjunction with
    the L{GougeConfigPrms} and L{GrainRandomBoxPacker} classes.
    Generates grains which consist of a particles packed into
    a circle/sphere. Spherical particle radii are generated
    independent of grain size.
    @see: getGrain
    """
    def __init__(
        self,
        minGRadius,
        maxGRadius,
        do2d=False,
        minPRadius=None,
        maxPRadius=None,
        maxInsertFails=5000
    ):
        """
        Initialise the generator.
        @type minGRadius: float
        @param minGRadius: minimum grain radius
        @type maxGRadius: float
        @param maxGRadius: maximum grain radius
        @type do2d: bool
        @param do2d: If True, grain particles all occur
          on the M{z=0} plane (ie particles packed into
          a circle).
        @type minPRadius: float
        @param minPRadius: The minimum radius of randomly
          generated particles.
        @type maxPRadius: float
        @param maxPRadius: The maximum radius of randomly
          generated particles.
        """
        if (maxPRadius == None):
            maxPRadius = maxGRadius*0.4
        if (minPRadius == None):
            minPRadius = maxPRadius*0.15

        self.do2d            = do2d
        self.maxInsertFails  = maxInsertFails
        self.minPRadius      = minPRadius
        self.maxPRadius      = maxPRadius
        RndGrainGen.__init__(
            self,
            minGRadius,
            maxGRadius,
            minPRadius,
            maxPRadius
        )

    def getGrain(self, p):
        """
        Uses L{RandomSpherePacker} to create a grain consisting
        of randomly sized particles packed into a sphere with
        centre C{p.getPosn()} and radius C{p.getRadius()}. The
        random particle radii will lie in the range
        C{(self.minPRadius,self.maxPRadius)}
        @type p: L{SimpleSphere}
        @param p: C{p.getPosn()} and C{p.getRadius()} determine
          the size and position of the returned grain.
        @rtype: L{Grain}
        @return: A I{spherical} grain of randomly sized particles.
        """
        minRadius = self.minPRadius
        maxRadius = self.maxPRadius
        packer = \
          RandomSpherePacker(
            minRadius=minRadius,
            maxRadius=maxRadius,
            cubicPackRadius=maxRadius,
            maxInsertFails=self.maxInsertFails,
            bSphere=BoundingSphere(p.getPosn(), p.getRadius()),
            do2d=self.do2d,
            tolerance = 0.1*minRadius
          )
        packer.generate()
        g = Grain()
        for p in packer.getParticleIterator():
            g.createParticle(p)
        getLogger().info(
          "Created grain of " + str(g.getNumParticles()) + " particles"
        )
        return g
