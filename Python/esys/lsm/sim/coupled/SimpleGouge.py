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
from __future__      import division
from esys.finley     import ReadMesh
from esys.linearPDEs import LinearPDE
from numarray        import identity
from esys.escript    import *
from LsmEscriptPy    import *

from esys.lsm.util   import Vec3

import re
import os
import logging

logger = logging.getLogger("esys.lsm.sim.coupled.SimpleGouge")
def getLogger():
    return logger

class NeighbourListPrms:
    def __init__(
        self,
        gridSpacing,
        updateDisplDist
    ):
        self.gridSpacing = gridSpacing
        self.updateDisplDist = updateDisplDist

class LsmGougePrms:
    def __init__(
        self,
        numWorkerProcesses,
        mpiDimList,
        lsmGeomFileName,
        meshFileName,
        timeStepSize,
        neighbourListPrms,
        bondPrms,
        frictionPrms,
        dampingPrms,
        bondedMeshPrms,
        chkptPrms,
        bBox = None
    ):
        self.numWorkerProcesses = numWorkerProcesses
        self.mpiDimList = mpiDimList
        self.lsmGeomFileName = lsmGeomFileName
        self.meshFileName = meshFileName
        self.timeStepSize = timeStepSize
        self.neighbourListPrms = neighbourListPrms
        self.bondPrms = bondPrms
        self.frictionPrms = frictionPrms
        self.dampingPrms = dampingPrms
        self.bondedMeshPrms = bondedMeshPrms
        self.chkptPrms = chkptPrms
        self.dim = None
        self.bBox = bBox

    def getNumWorkerProcesses(self):
        return self.numWorkerProcesses

    def getMpiDimList(self):
        return self.mpiDimList

    def getMeshFileName(self):
        return self.meshFileName

    def getMeshName(self):
        return self.bondedMeshPrms.getMeshName()

    def getDim(self):
        if (self.dim == None):
            regex = re.compile("^Dimension\s+(\d)D")
            f = file(self.lsmGeomFileName, "r")
            line = f.readline()
            match = regex.match(line)
            while (len(line) > 0) and (match == None):
                line = f.readline()
                match = regex.match(line)
            f.close()
            if (match != None):
                self.dim = int(match.group(1))
            else:
                raise Exception(
                    "Could not locate dimension info in file " \
                    + \
                    str(self.lsmGeomFileName)
                )
        return self.dim

    def is2d(self):
        return (self.getDim() == 2)

    def initialise(self, lsm):
        lsm.initVerletModel(
            self.getParticleType(),
            self.neighbourListPrms.gridSpacing,
            self.neighbourListPrms.updateDisplDist,
        )
        lsm.setTimeStepSize(self.timeStepSize)
        if (self.bBox != None):
            lsm.setSpatialDomain(self.bBox.getMinPt(), self.bBox.getMaxPt())
        lsm.readGeometry(self.lsmGeomFileName)
        if (self.is2d()):
            lsm.force2dComputations()
        lsm.createInteractionGroup(self.bondPrms)
        lsm.createInteractionGroup(self.frictionPrms)
        lsm.createExclusion(self.bondPrms.getName(), self.frictionPrms.getName())
        lsm.createDamping(self.dampingPrms)
        if (self.is2d()):
            lsm.readMesh2D(self.meshFileName, self.getMeshName(), 20)
        else:
            lsm.readMesh(self.meshFileName, self.getMeshName())
        lsm.createInteractionGroup(self.bondedMeshPrms)
        lsm.createCheckPointer(self.chkptPrms)

class LsmSimpleGouge(LsmMpiEscript):
    def __init__(self, gougePrms):
        self.gougePrms = gougePrms
        LsmMpiEscript.__init__(
            self,
            self.gougePrms.getNumWorkerProcesses(),
            self.gougePrms.getMpiDimList()
        )
        self.initialise()

    def initialise(self):
        self.gougePrms.initialise(self)

    def is2d(self):
        return self.gougePrms.is2d()

    def getDim(self):
        return self.gougePrms.getDim()

    def getMeshFileName(self):
        return self.gougePrms.getMeshFileName()

    def getMeshName(self):
        return self.gougePrms.getMeshName()

    def visitRefStressPairs(self, visitor):
        if (self.is2d):
            LsmMpiEscript.visitRefStressPairs2d(self, self.getMeshName(), visitor)
        else:
            LsmMpiEscript.visitRefStressPairs(self, self.getMeshName(), visitor)

    def visitNodeRefs(self, visitor):
        if (self.is2d):
            LsmMpiEscript.visitNodeRefs2d(self, self.getMeshName(), visitor)
        else:
            LsmMpiEscript.visitNodeRefs(self, self.getMeshName(), visitor)

    def moveNodeBy(self, nodeRef, displacement):
        self.moveSingleMeshNodeBy(self.getMeshName(), nodeRef, displacement)

class EscriptGougePrms:
    def __init__(
        self,
        lameLambda, # Lame coeff
        lameMu,     # Lame coeff
        density      = 1.0,
        viscosity    = 0.0,
        pressure     = 0.001,
        shearForce   = 0.005,
        maxTime      = 1000,
        timeStepSize = 0.05,
        saveDxIncr   = 10,
        outputDir    = "."
    ):
        self.lameLambda   = lameLambda
        self.lameMu       = lameMu
        self.density      = density
        self.viscosity    = viscosity
        self.pressure     = pressure
        self.shearForce   = shearForce
        self.maxTime      = maxTime
        self.timeStepSize = timeStepSize
        self.outputDir    = outputDir
        self.saveDxIncr   = saveDxIncr

class NodeRefVisitor:
    def __init__(self, data, lsm):
        self.data = data
        self.lsm  = lsm

    def visitNodeRef(self, nodeRef):
        displArray = numarray.array([0.0]*self.lsm.getDim())
        self.data.getRefValue(nodeRef, displArray)
        displVec = Vec3()
        for i in range(0, self.lsm.getDim()):
            displVec[i] = displArray[i]
        getLogger().debug("Moving node " + str(nodeRef) + " by " + str(displVec))
        self.lsm.moveNodeBy(nodeRef, displVec)

class RefStressVisitor:
    def __init__(self, data, dim):
        self.data = data
        self.dim = dim

    def visitRefStressPair(self, ref, stressVec):
        getLogger().debug("Setting stress on boundary element " + str(ref) + " to " + str(stressVec))
        self.data.setRefValue(
            ref,
            numarray.array(stressVec.toList()[0:self.dim])
        )

class EscriptSimpleGouge:
    def __init__(self, lsm, prms=None):
        self.lsm = lsm
        self.prms = prms
        if (self.prms == None):
            self.initialisePrmsFromLsm()

    def initialisePrmsFromLsm(self):
        pass

    def getDim(self):
        return self.lsm.getDim()

    def getMeshFileName(self):
        return self.lsm.getMeshFileName()

    def updateImpactStresses(self, stressData):
        if (self.lsm != None):
            self.lsm.visitRefStressPairs(
                RefStressVisitor(data=stressData, dim=self.lsm.getDim())
            )

    def updateLsmMeshPosn(self, displacementData):
        if (self.lsm != None):
            self.lsm.visitNodeRefs(NodeRefVisitor(displacementData, self.lsm))

    def run(self):
        # material parameter
        lam = self.prms.lameLambda
        mu  = self.prms.lameMu
        rho = self.prms.density
        nu  = self.prms.viscosity
        
        # set up boundary conditions
        pres  = self.prms.pressure
        shear = self.prms.shearForce
        
        getLogger().info("Reading mesh from " + self.getMeshFileName())
        domain=ReadMesh(self.getMeshFileName())
        impact_forces=escript.Vector(0,FunctionOnBoundary(domain))
        impact_forces.expand()
        getLogger().info("Shape = " + str(impact_forces.getShape()))
        x=FunctionOnBoundary(domain).getX()
        getLogger().info("Initialising pressure and shearing boundary conditions...")
        snapDist = 0.001
        external_forces = \
            (abs(x[1]-sup(x[1]))-snapDist).whereNegative()*[-shear,-pres] \
            +                                                          \
            (abs(x[1]-inf(x[1]))-snapDist).whereNegative()*[shear,pres]

        getLogger().info("Setting up PDE...")
        mypde = LinearPDE(domain)
        mypde.setLumpingOn()
        mypde.setValue(D=rho*identity(mypde.getDim()))
        
        getLogger().info("Initialising solution at t=0...")
        u      = Vector(0,ContinuousFunction(domain))
        u_last = Vector(0,ContinuousFunction(domain))
        v      = Vector(0,ContinuousFunction(domain))
        v_last = Vector(0,ContinuousFunction(domain))
        a      = Vector(0,ContinuousFunction(domain))
        a_last = Vector(0,ContinuousFunction(domain))
        
        # initialise iteration prms
        tend = self.prms.maxTime
        dt   = self.prms.timeStepSize
        # dt=1./5.*sqrt(rho/(lam+2*mu))*Lsup(domain.getSize())
        getLogger().info("time step size = " + str(dt))
        n=0
        t=0
        
        getLogger().info("Beginning iteration...")
        while (t < tend):
            getLogger().info("Running LSM time step...");
            self.lsm.runTimeStep()
            
            getLogger().info("Updating impact forces from LSM...");
            self.updateImpactStresses(impact_forces)
            getLogger().info(
              "(inf(impactForces), sup(impact_forces)) = (" + \
              str(inf(impact_forces)) + ", " + str(sup(impact_forces)) + ")"
            )
            
            
            # ... update FEM ...
            getLogger().info("Initialising PDE coefficients...")
            g=grad(u)
            stress=(lam*trace(g))*identity(mypde.getDim())+mu*(g+transpose(g))
            mypde.setValue(
                X = -(1.0/(1.0+(nu*dt/2.0)))*stress,
                Y = - nu*v_last-(nu*dt/2.0)*a_last,
                y=external_forces+impact_forces
            )
            getLogger().info("Solving PDE...")
            a=mypde.getSolution()
            
            getLogger().info("Updating displacements...")
#            u_new=2*u-u_last+dt**2*a
#            u_last=u
#            u=u_new

            v_new = v_last + (dt/2.0)*(a + a_last)
            v_last = v
            v = v_new

            u_new = u_last + dt*v + ((dt**2)/2.0)*a
            u_last = u
            u = u_new
            
            a_last = a

            getLogger().info("Updating LSM mesh node positions...")
            displacement = u - u_last
            self.updateLsmMeshPosn(displacement)
            
            t+=dt
            n+=1
            getLogger().info(str(n) + "-th time step, t=" + str(t))
            getLogger().info("a=" + str(inf(a)) + ", " + str(sup(a)))
            getLogger().info("u=" + str(inf(u)) + ", " + str(sup(u)))
            getLogger().info("inf(u-u_last) = " + str(inf(displacement)))
            getLogger().info("sup(u-u_last) = " + str(sup(displacement)))
            
            # ... save current acceleration and displacement
            if ((self.prms.saveDxIncr > 0) and ((n%self.prms.saveDxIncr) == 0)):
                u.saveDX(os.path.join(self.prms.outputDir, "displ.{0:d}.dx".format(n//self.prms.saveDxIncr)))

class SimpleGouge:
    def __init__(
        self,
        lsmGougePrms,
        escriptPrms
    ):
        self.lsmPrms = lsmGougePrms
        self.escriptPrms = escriptPrms

    def run(self):
        lsm = LsmSimpleGouge(self.lsmPrms)
        escriptSim = EscriptSimpleGouge(lsm, self.escriptPrms)
        escriptSim.run()

