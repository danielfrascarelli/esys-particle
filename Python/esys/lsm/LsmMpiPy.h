/////////////////////////////////////////////////////////////
//                                                         //
// Copyright (c) 2003-2017 by The University of Queensland //
// Centre for Geoscience Computing                         //
// http://earth.uq.edu.au/centre-geoscience-computing      //
//                                                         //
// Primary Business: Brisbane, Queensland, Australia       //
// Licensed under the Open Software License version 3.0    //
// http://www.apache.org/licenses/LICENSE-2.0              //
//                                                         //
/////////////////////////////////////////////////////////////

#ifndef ESYS_LSM_LSMMPIPY_H
#define ESYS_LSM_LSMMPIPY_H

// --- boost includes ---
#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>

// --- STL includes ---
#include <string>
#include <vector>
#include <map>

// --- Project includes ---
#include "Parallel/LatticeMaster.h"
#include "Python/esys/lsm/InteractionGroupPy.h"
#include "Python/esys/lsm/BondInteractionGroupPy.h"

using std::string;

namespace esys
{
  namespace lsm
  {
    void exportLsm();

    class ParticleScalarFieldSaverPrmsPy;
    class ParticleVectorFieldSaverPrmsPy;
    class TaggedParticleScalarFieldSaverPrmsPy;
    class TaggedParticleVectorFieldSaverPrmsPy;
    class InteractionScalarFieldSaverPrmsPy;
    class CheckedInteractionScalarFieldSaverPrmsPy;
    class TaggedInteractionScalarFieldSaverPrmsPy;
    class InteractionVectorFieldSaverPrmsPy;
    class CheckedInteractionVectorFieldSaverPrmsPy;
    class WallVectorFieldSaverPrmsPy;
    class TriangleScalarFieldSaverPrmsPy;
    class TriangleVectorFieldSaverPrmsPy;
    class NRotBondPrmsPy;
    class CappedNRotBondPrmsPy;
    class NRotShortBondPrmsPy;
    class NRotElasticPrmsPy;
    class HertzianElasticPrmsPy;
    class HertzianViscoElasticFrictionPrmsPy;
    class HertzianViscoElasticPrmsPy;
    class HertzMindlinPrmsPy;
    class HertzMindlinViscoPrmsPy;
    class LinearDashpotPrmsPy;
    class NRotFrictionPrmsPy;
    class SpringDashpotFrictionPrmsPy;
    class RotBondPrmsPy;
    class RotThermBondPrmsPy;
    class BrittleBeamPrmsPy;
    class BrittleBeamSCPrmsPy;
    class BrittleBeamDZCPrmsPy;
    class FrictionPrmsPy;
    class RotFrictionPrmsPy;
    class RotThermFrictionPrmsPy;
    class RotElasticPrmsPy;
    class RotThermElasticPrmsPy;
    class VWFrictionPrmsPy;
    class DampingPrmsPy;
    class LocalDampingPrmsPy;
    class RotLocalDampingPrmsPy;
    class ABCDampingPrmsPy;
    class NRotBondedWallPrmsPy;
    class NRotSoftBondedWallPrmsPy;
    class NRotElasticWallPrmsPy;
    class NRotElasticSphereBodyPrmsPy;
    class CheckPointPrmsPy;
    class RestartCheckPointPrmsPy;
    class NRotElasticTriMeshPrmsPy;
    class NRotElasticMesh2DPrmsPy;
    class NRotElasticLinMeshPrmsPy;
    class NRotBondedTriMeshPrmsPy;
    class NRotBondedLinMeshPrmsPy;
    class RunnablePy;
    class Vec3Py;
    class BoundingBoxPy;
    class GravityPrmsPy;
    class BuoyancyPrmsPy;
    class MaxTriggerPrmsPy;
    class FluidScalarFieldSaverPrmsPy; //fluid contents
    class FluidVectorFieldSaverPrmsPy; //fluid contents
    class FluidInteractionScalarFieldSaverPrmsPy; //fluid contents
    class FluidInteractionVectorFieldSaverPrmsPy; //fluid contents


    void checkMpiDimensions(int numProcesses, const std::vector<int> &mpiDimVector);

    void checkMpiDimensionsPy(int numProcesses, const boost::python::list &mpiDimList);

    void checkParticleType(const std::string &particleType);

    /*!
      \brief Wrapper to make LatticeMaster methods available in Python
    */
    class LsmMpiPy
    {
    public:
      typedef CLatticeMaster::ParticleIdPair ParticleIdPair;
      typedef CLatticeMaster::ParticleIdPairVector ParticleIdPairVector;
      typedef CLatticeMaster::MeshNodeDataVector MeshNodeDataVector;
      typedef CLatticeMaster::MeshTriDataVector  MeshTriDataVector;
      typedef CLatticeMaster::TriMeshDataPair    TriMeshDataPair;



      LsmMpiPy(
        int numWorkerProcesses,
        const boost::python::list &mpiDimList
      );

      virtual ~LsmMpiPy();

      int getNumWorkerProcesses() const;

      void initVerletModel(
        const std::string &particleType,
        double gridSpacing,
        double verletDist
      );

      double getTimeStepSize() const;

      void setTimeStepSize(double dt);

      void setTimingFileName(const std::string &fileNamePrefix);
      void setSlaveTimingFileName(const std::string &fileNamePrefix);

      std::string getParticleType() const;
      std::string getLsmVersion() const;
      void readGeometry(const std::string &fileName);

      int getNumParticles();

      int getTimeStep() const;

      void createParticles(boost::python::object &iterable);

      void createParticle(boost::python::object &particle);

      void createConnections(boost::python::object &iterable);

      /****fluid contents: begin****/
      void createFluidInteraction(double,double,double,double,double,double,double,const Vec3Py&,const Vec3Py&,double);
      void createFluidInteractionVec3(const Vec3Py&,double,double,double,double,double,double,const Vec3Py&,const Vec3Py&,double);
      void createFluidScalarFieldSaver(const FluidScalarFieldSaverPrmsPy &prms);
      void createFluidVectorFieldSaver(const FluidVectorFieldSaverPrmsPy &prms);
      void createFluidInteractionScalarFieldSaver(const FluidInteractionScalarFieldSaverPrmsPy &prms);
      void createFluidInteractionVectorFieldSaver(const FluidInteractionVectorFieldSaverPrmsPy &prms);
      /****fluid contents: end****/

      // --- interaction creation functions ---
      void createNRotElasticInteractGrp(const NRotElasticPrmsPy &prms);
      void createNRotBondInteractGrp(const NRotBondPrmsPy &bondPrms);
      void createCappedNRotBondInteractGrp(const CappedNRotBondPrmsPy &bondPrms);
      void createNRotShortBondInteractGrp(const NRotShortBondPrmsPy &bondPrms);
      void createNRotFrictionInteractGrp(const NRotFrictionPrmsPy &prms);
      void createSpringDashpotFrictionInteractGrp(const SpringDashpotFrictionPrmsPy &prms);
      void createRotBondInteractGrp(const RotBondPrmsPy &bondPrms);
      BondInteractionGroupPy createRotThermBondInteractGrp(const RotThermBondPrmsPy &bondPrms);
      void createBrittleBeamInteractGrp(const BrittleBeamPrmsPy &bondPrms);
      void createBrittleBeamSCInteractGrp(const BrittleBeamSCPrmsPy &bondPrms);
      void createBrittleBeamDZCInteractGrp(const BrittleBeamDZCPrmsPy &bondPrms);
      void createFrictionInteractGrp(const FrictionPrmsPy &prms);
      void createRotFrictionInteractGrp(const RotFrictionPrmsPy &prms);
      void createRotThermFrictionInteractGrp(const RotThermFrictionPrmsPy &prms);
      void createRotElasticInteractGrp(const RotElasticPrmsPy &prms);
      void createRotThermElasticInteractGrp(const RotThermElasticPrmsPy &prms);
      void createDamping(const DampingPrmsPy &prms);
      void createLocalDamping(const LocalDampingPrmsPy &prms);
      void createRotLocalDamping(const RotLocalDampingPrmsPy &prms);
      void createABCDamping(const ABCDampingPrmsPy &prms);
      void createGravity(const GravityPrmsPy&);
      void createBuoyancy(const BuoyancyPrmsPy&);
      void createVWFrictionIG(const VWFrictionPrmsPy&);
      void createHertzianElasticIG(const HertzianElasticPrmsPy &prms);
      void createHertzianViscoElasticFrictionIG(const HertzianViscoElasticFrictionPrmsPy &prms);
      void createHertzianViscoElasticIG(const HertzianViscoElasticPrmsPy &prms);
      void createHertzMindlinIG(const HertzMindlinPrmsPy &prms);
      void createHertzMindlinViscoIG(const HertzMindlinViscoPrmsPy &prms);
      void createLinearDashpotIG(const LinearDashpotPrmsPy &prms);


      // --- remove interactions ---
      void removeInteractionGrp(const std::string&);

      // --- tagged interaction creation functions ---
      void createRotFrictionInteractGrpTag(const RotFrictionPrmsPy &prms,int,int,int,int);
      void createFrictionInteractGrpTag(const FrictionPrmsPy &prms,int,int,int,int);
      void createNRotFrictionInteractGrpTag(const NRotFrictionPrmsPy &prms,int,int,int,int);
      void createSpringDashpotFrictionInteractGrpTag(const SpringDashpotFrictionPrmsPy &prms,int,int,int,int);
      void createLinearDashpotInteractGrpTag(const LinearDashpotPrmsPy &prms,int,int,int,int);
      void createRotElasticInteractGrpTag(const RotElasticPrmsPy &prms,int,int,int,int);
      void createElasticInteractGrpTag(const NRotElasticPrmsPy &prms,int,int,int,int);


      void createExclusion(
        const std::string &interactionName1,
        const std::string &interactionName2
      );


      // --- particle property setting functions ---
      void setParticleVel(int,const Vec3Py&);
      void setParticleAngVel(int,const Vec3Py&);
      void setParticleDensity(int,int,double);
      void resetParticleOrientation(int,int);
      void setTaggedParticleVel(int,const Vec3Py&);
      void setVelocityOfWall(const std::string&,const Vec3Py&);
      void tagParticleNearestTo(int,int,const Vec3Py&);
      void setParticleNonDynamic(int);
      void setParticleNonRot(int);
      void setParticleNonTrans(int);
      
      // --- modifying interaction parameters ---
      void setInteractionParameter(const std::string&,const std::string&,double);
      
      // ---- checkpointing -----------
      void createCheckPointer(const RestartCheckPointPrmsPy &prms);
      void createCheckPointerThroughMaster(const RestartCheckPointPrmsPy &prms); // write through master
      void createSnapShots(const CheckPointPrmsPy &prms);
      void loadCheckPoint(const std::string&);

      // ------------------------------
      int getNumTimeSteps() const;

      void setNumTimeSteps(int numTimeSteps);

      // --- Mesh functions ---

      void readMeshWithTag(const std::string &fileName, const std::string &meshName, int tag);
      void readMesh(const std::string &fileName, const std::string &meshName);
      void createTriMesh(
        const std::string &meshName,
        const boost::python::object &nodeSequence,
        const boost::python::object &triSequence
      );
      void translateMesh(const std::string&,const Vec3Py&);
      void rotateMesh(const std::string&,const Vec3Py&,const Vec3Py&,double);
      void readMesh2D(const std::string &fileName, const std::string &meshName, int tag);

      void createNRotElasticTriMeshInteractGrp(const NRotElasticTriMeshPrmsPy &prms);

      void createNRotBondedTriMeshInteractGrp(const NRotBondedTriMeshPrmsPy &prms);
      void createNRotElasticMesh2DInteractGrp(const NRotElasticMesh2DPrmsPy &prms);
      void createNRotElasticLinMeshInteractGrp(const NRotElasticLinMeshPrmsPy &prms);
      void createNRotBondedLinMeshInteractGrp(const NRotBondedLinMeshPrmsPy &prms);

      void moveSingleMeshNodeBy(const std::string& meshname, int id, const Vec3Py& d);

      void addPreTimeStepRunnable(RunnablePy &runnable);

      void addPostTimeStepRunnable(RunnablePy &runnable);

      void force2dComputations(bool do2d);

      void setBBoxSpatialDomain(const BoundingBoxPy &domain);

      void setBBoxSpatialDomainWithCirc(
        const BoundingBoxPy &domain,
        const boost::python::list &circDimList
      );

      void setSpatialDomain(const Vec3Py &minPt, const Vec3Py &maxPt);

      void getInitMinMaxPt(Vec3Py &initMinPt, Vec3Py &initMaxPt);

      void getCurMinMaxPt(Vec3Py &curMinPt, Vec3Py &curMaxPt);

      int findClosestParticle(const Vec3Py &pt);

      Vec3Py getParticlePosn(int particleId);

      void changeRadiusBy(int tag, const double &deltaR);
      // --- move particles ---
      void moveTaggedParticlesTo(int tag, const Vec3Py &pt);
      void moveTaggedParticlesBy(int tag, const Vec3Py &displacement);
      void moveSingleParticleTo(int particleId, const Vec3Py &pt);

      // --- wall related functions ---
      void createWall(const string &name, const Vec3Py &posn, const Vec3Py &normal);
      void createSphereBody(const string &name, const Vec3Py &posn, const double &radius);
      void createNRotBondedWall(const NRotBondedWallPrmsPy &prms);
      void createNRotElasticWall(const NRotElasticWallPrmsPy &prms);
      void createNRotElasticSphereBody(const NRotElasticSphereBodyPrmsPy &prms);
      void createNRotSoftBondedWall(const NRotSoftBondedWallPrmsPy &prms);
      void createNRotElasticWallTagged(const NRotElasticWallPrmsPy &prms, int tag, int mask);
      void moveWallBy(const string&, const Vec3Py &disp);
      void moveSphereBodyBy(const string&, const Vec3Py &disp);
      void setWallNormal(const string&, const Vec3Py &wn);
      void applyForceToWall(const string&, const Vec3Py&);
      Vec3Py getWallPosition(const std::string&);
      Vec3Py getWallForce(const std::string&);
      Vec3Py getSphereBodyPosition(const std::string&);
      Vec3Py getSphereBodyForce(const std::string&);


      void runTimeStep();
      void run();

      // Exit the simulation after running a series of single steps
      // of the time-integration method.
      void exit();

      // --- console related functions ---
      void SetVerbosityPy(bool);
      void SetVerbosityLevelPy(int);
      void SetConsoleFilenamePy(const std::string&);
      void SetConsoleBufferedPy(unsigned int);

      // --- field saving functions ---
      void createParticleScalarFieldSaver(
        const ParticleScalarFieldSaverPrmsPy &prms
      );

      void createParticleVectorFieldSaver(const ParticleVectorFieldSaverPrmsPy&);
      void createInteractionScalarFieldSaver(const InteractionScalarFieldSaverPrmsPy &prms);
      void createCheckedInteractionScalarFieldSaver(const CheckedInteractionScalarFieldSaverPrmsPy &prms);
      void createInteractionVectorFieldSaver(const InteractionVectorFieldSaverPrmsPy&);
      void createCheckedInteractionVectorFieldSaver(const CheckedInteractionVectorFieldSaverPrmsPy &prms);

      void createTaggedParticleScalarFieldSaver(const TaggedParticleScalarFieldSaverPrmsPy&);
      void createTaggedParticleVectorFieldSaver(const TaggedParticleVectorFieldSaverPrmsPy&);
      void createTaggedInteractionScalarFieldSaver(const TaggedInteractionScalarFieldSaverPrmsPy&);

      void addTaggedScalarParticleDistributionSaver(const string&,const string&,const string&,int,int,int,int,int,int,double,double,int);
      void addVectorTriangleSaveField(const TriangleVectorFieldSaverPrmsPy&);
      void addScalarTriangleSaveField(const TriangleScalarFieldSaverPrmsPy&);
      void addVectorWallField(const WallVectorFieldSaverPrmsPy &prms);

      // --- fields with trigger ---
      void createParticleVectorFieldSaverWithTrigger(const MaxTriggerPrmsPy&,const ParticleVectorFieldSaverPrmsPy&);
      void createTaggedParticleVectorFieldSaverWithTrigger(const MaxTriggerPrmsPy&,const TaggedParticleVectorFieldSaverPrmsPy&);

      void visitNodeRefs2d(const std::string &meshName, boost::python::object pyObject);

      void visitRefStressPairs2d(const std::string &meshName, boost::python::object pyObject);

      void visitRefForcePairs(const std::string &meshName, boost::python::object pyObject);

      void visitParticlesWithId(
        const boost::python::list &idList,
        boost::python::object &pyObject
      );

      void visitParticles(
        boost::python::object &pyObject
      );

      boost::python::list getParticleList();

      boost::python::list getParticleWithIdList(
        const boost::python::list &idList
      );

      boost::python::list getParticlesInBBox(
        const BoundingBoxPy &bbox
      );

      void createBonds(
        const std::string          &groupName,
        const ParticleIdPairVector &idPairVector
      );

      void updateInteractions();

      ParticleIdPairVector getBondGroupIdPairs(
        const std::string &groupName
      );
      void setVerbosityPy(int);


    protected:
      typedef std::map<std::string, std::string> InteractionNameTypeMap;

      InteractionNameTypeMap &getNameTypeMap();
      const InteractionNameTypeMap &getNameTypeMap() const;

      const CLatticeMaster &getLatticeMaster() const;

      CLatticeMaster &getLatticeMaster();

    private:
      class Impl;
      typedef boost::shared_ptr<Impl> ImplPtr;
      ImplPtr m_implPtr;
    };


    void setVerbosityPy(bool);
    void setVerbosityLevelPy(int);
  }
}
#endif

