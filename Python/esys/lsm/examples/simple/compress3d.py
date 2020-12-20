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
Example simulation of compression of an elastically bonded block of
non-rotational spheres.
"""
from esys.lsm      import *
from esys.lsm.util import Vec3, BoundingBox, InstallInfo

class Loading(Runnable):
    """
    Loading mechanism which moves the walls.
    """
    def __init__(self,lsm):
        """
        Initialise loading object with a reference to an L{LsmMpi}
        object.
        @type lsm: L{esys.lsm.LsmMpi<esys.lsm.LsmPy.LsmMpi>}
        @param lsm: Model object which will move the model wall.
        """
        Runnable.__init__(self)
        self.theLSM=lsm

    def run(self):

        """
        Moves walls in compressive fashion.
        """
        self.theLSM.moveWallBy("lowerWall",Vec3(0.0,0.0001,0.0))
        self.theLSM.moveWallBy("upperWall",Vec3(0.0,-0.0001,0.0))

def runSimulation():
    """
    Runs a compression simulation on an elastic block.
    Outputs wall-forces and wall-positions to file.
    """
    setVerbosity(True)
    mySim=LsmMpi(2,[0,0,0])
    mySim.initVerletModel("NRotSphere", 2.5, 0.5)
    mySim.setTimeStepSize(0.02)
    mySim.setSpatialDomain(BoundingBox(Vec3(-5.0,0.0,-5.0),Vec3(15.0,10.0,15.0)))
    mySim.readGeometry(
      InstallInfo.getDataFilePath("cube10r0.2.geo")
    );
    # setup interactions
    bip=NRotBondPrms(1,"bonded",1.0,1.05)
    fip=NRotFrictionPrms("friction",1.0,0.6,1.0)
    mySim.createInteractionGroup(bip)
    mySim.createInteractionGroup(fip)
    mySim.createExclusion("bonded","friction")
    # wall parameters
    mySim.createWall("lowerWall",Vec3(0.0,0.0,0.0),Vec3(0.0,1.0,0.0))
    mySim.createWall("upperWall",Vec3(0.0,10.0,0.0),Vec3(0.0,-1.0,0.0))
    wp1=NRotElasticWallPrms("upperWallInteraction","upperWall",1.0)
    wp2=NRotElasticWallPrms("lowerWallInteraction","lowerWall",1.0)
    # setup savers
    mySim.createFieldSaver(
      WallVectorFieldSaverPrms(
        fileName="wf2.dat",
        fieldName="Force",
        wallName=["lowerWall","upperWall"],
        fileFormat="RAW_SERIES",
        beginTimeStep=0,
        endTimeStep=10,
        timeStepIncr=1
      )
    )
    mySim.createFieldSaver(
      WallVectorFieldSaverPrms(
        fileName="wp2.dat",
        fieldName="Position",
        wallName=["lowerWall","upperWall"],
        fileFormat="RAW_SERIES",
        beginTimeStep=0,
        endTimeStep=10,
        timeStepIncr=1
      )
    )
    # create walls
    mySim.createInteractionGroup(wp1)
    mySim.createInteractionGroup(wp2)
    # add loading function
    lf=Loading(mySim)
    mySim.addPreTimeStepRunnable(lf)
    mySim.setNumTimeSteps(10)
    mySim.run()

if (__name__=="__main__"):
    runSimulation()
