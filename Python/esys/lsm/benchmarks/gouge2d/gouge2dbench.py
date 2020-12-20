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
Defines L{runSimulation} function which executes gouge shear simulation.
"""
from __future__    import division, print_function
from time          import *
import os
import sys
import popen2

from esys.lsm      import *
from esys.lsm.util import Vec3, BoundingBox, InstallInfo, pathSearcher

class Loading(Runnable):
    """
    Objects of this class provide the loading mechanism for a compression
    simulation.
    """
    def __init__(self,lsm,dt,v):
        Runnable.__init__(self)
        self.theLSM=lsm
        self.dx=v*dt

    def run(self):
        """
        Moves upper and lower walls by an increment parallel to the gouge
        axis and applies a constant force perpendicular to it.
        """
        t=self.theLSM.getTimeStep()
        if(t<1000):
            dx_cur=self.dx*(t/1000.0)
        else:
            dx_cur=self.dx
        self.theLSM.moveWallBy("lowerWall",Vec3(dx_cur,0.0,0.0))
        self.theLSM.moveWallBy("upperWall",Vec3(-1.0*dx_cur,0.0,0.0))
        self.theLSM.applyForceToWall("upperWallInteraction",Vec3(0.0,-0.25,0.0));
        self.theLSM.applyForceToWall("lowerWallInteraction",Vec3(0.0,0.25,0.0));
        
def runSimulation():
    """
    Initialises elastic block model and runs the compression simulation.
    """
    dt=0.02
    v=0.002
    nt=10000
    ncpu_x=int(sys.argv[1])
    ncpu_y=int(sys.argv[2])
##    setVerbosity(True)
    mySim=LsmMpi(ncpu_x*ncpu_y,[ncpu_x,ncpu_y,1])
    mySim.initVerletModel("RotSphere", 2.5, 0.5)
    mySim.setTimeStepSize(dt)
    mySim.force2dComputations(True)
    # read geometry
    mySim.readGeometry("bench_gouge.geo")
    # setup interactions
    bip=RotBondPrms(0,"bonded",0.5,0.15,0.04,0.017,0.025,0.125,0.0125,0.0125)
    fip=RotFrictionPrms("friction", 1.0, 0.6, 0.6, 1.0)
    dip=DampingPrms("Damping","damping1",0.01,50)
    rdip=DampingPrms("RotDamping","damping2",0.01,50)
    mySim.createInteractionGroup(bip)
    mySim.createInteractionGroup(fip)
    mySim.createExclusion("bonded","friction")
    mySim.createInteractionGroup(dip)
    mySim.createInteractionGroup(rdip)
    # create walls
    mySim.createWall("lowerWall",Vec3(0.0,0.0,0.0),Vec3(0.0,1.0,0.0))
    mySim.createWall("upperWall",Vec3(0.0,40.0,0.0),Vec3(0.0,-1.0,0.0))
    wp1=NRotBondedWallPrms("upperWallInteraction","upperWall",1.0,4)
    wp2=NRotBondedWallPrms("lowerWallInteraction","lowerWall",1.0,3);
    mySim.createInteractionGroup(wp1)
    mySim.createInteractionGroup(wp2)
    # setup savers
    mySim.createFieldSaver(
      WallVectorFieldSaverPrms(
        fileName="wf.dat",
        fieldName="Force",
        wallName=["lowerWall","upperWall"],
        fileFormat="RAW_SERIES",
        beginTimeStep=0,
        endTimeStep=nt,
        timeStepIncr=10
      )
    )
    mySim.createFieldSaver(
      WallVectorFieldSaverPrms(
        fileName="wp.dat",
        fieldName="Position",
        wallName=["lowerWall","upperWall"],
        fileFormat="RAW_SERIES",
        beginTimeStep=0,
        endTimeStep=nt,
        timeStepIncr=10
      )
    )
    ek_prm=ParticleScalarFieldSaverPrms("e_kin","ekin.dat","SUM",0,nt,10)
    mySim.createFieldSaver(ek_prm)
    # add loading function
    lf=Loading(mySim,dt,v)
    mySim.addPreTimeStepRunnable(lf)
    mySim.setNumTimeSteps(nt)
    nparts=mySim.getNumParticles()
    print("Particles in the model: ", nparts)
    start_time=time()
    mySim.run()
    stop_time=time()
    run_time=stop_time-start_time
    # calculate relative performance
    perf=(nparts*nt)/run_time
    # print results
    print("runtime     : ", run_time, " seconds") 
    print("performance : ", perf, " particles*timesteps/second")

if (__name__ == "__main__"):
    runSimulation()
