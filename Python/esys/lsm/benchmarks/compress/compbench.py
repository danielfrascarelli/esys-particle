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
Defines L{runSimulation} function which executes compression simulation.
"""
from __future__    import print_function

from esys.lsm      import *
from esys.lsm.util import Vec3, BoundingBox, InstallInfo
from time          import *

class Loading(Runnable):
    """
    Objects of this class provide the loading mechanism for a compression
    simulation.
    """
    def __init__(self,lsm):
        Runnable.__init__(self)
        self.theLSM=lsm

    def run(self):
        """
	      Moves upper and lower walls by an increment.
	      """
        self.theLSM.moveWallBy("lowerWall",Vec3(0.0,0.00005,0.0))
        self.theLSM.moveWallBy("upperWall",Vec3(0.0,-0.00005,0.0))

def runSimulation():
    """
    Initialises elastic block model and runs the compression simulation.
    """
    #setVerbosity(True)
    mySim=LsmMpi(1,[0,0,0])
    mySim.initVerletModel("RotSphere", 2.5, 0.5)
    mySim.setTimeStepSize(0.01)
    mySim.setSpatialDomain(BoundingBox(Vec3(-5.0,-5.0,-5.0),Vec3(15.0,25.0,15.0)))
    mySim.readGeometry(
      InstallInfo.getDataFilePath("bench_block_10x20x10_r0.15.geo")
    );
    # setup interactions 
    bip=RotBondPrms(0,"bonded",0.5,0.15,0.04,0.017,0.0025,0.0125,0.00125,0.00125)
    fip=RotFrictionPrms( "friction", 1.0, 0.6, 0.6, 1.0)
    dip=DampingPrms("Damping","damping1",0.01,50)
    rdip=DampingPrms("RotDamping","damping2",0.01,50)
    mySim.createInteractionGroup(bip)
    mySim.createInteractionGroup(fip)
    mySim.createExclusion("bonded","friction")
    mySim.createInteractionGroup(dip)
    mySim.createInteractionGroup(rdip)
    # create walls
    mySim.createWall("lowerWall",Vec3(0.0,0.0,0.0),Vec3(0.0,1.0,0.0))
    mySim.createWall("upperWall",Vec3(0.0,20.0,0.0),Vec3(0.0,-1.0,0.0))
    wp1=NRotElasticWallPrms("upperWallInteraction","upperWall",1.0)
    wp2=NRotElasticWallPrms("lowerWallInteraction","lowerWall",1.0);
    mySim.createInteractionGroup(wp1)
    mySim.createInteractionGroup(wp2)
    # setup savers
    nb_prm=InteractionScalarFieldSaverPrms("bonded","count","nbonds","SUM",0,12000,4)
    mySim.createFieldSaver(nb_prm)
    mySim.createFieldSaver(
      WallVectorFieldSaverPrms(
        fileName="wf.dat",
        fieldName="Force",
        wallName=["lowerWall","upperWall"],
        fileFormat="RAW_SERIES",
        beginTimeStep=0,
        endTimeStep=12000,
        timeStepIncr=4
      )
    )
    # add loading function
    lf=Loading(mySim)
    mySim.addPreTimeStepRunnable(lf)
    mySim.setNumTimeSteps(12000)
    start_time=time()
    mySim.run()
    stop_time=time()
    print("runtime: ", stop_time-start_time, " seconds") 

if (__name__ == "__main__"):
    runSimulation()

