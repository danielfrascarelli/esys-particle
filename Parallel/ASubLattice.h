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

#ifndef __ASUBLATTICE_H
#define __ASUBLATTICE_H

// -- project includes --
#include "Model/EWallInteractionGroup.h"
#include "Model/BWallInteractionGroup.h"
#include "Model/SoftBWallInteractionGroup.h"
#include "Model/ViscWallIG.h"
#include "Parallel/CheckPointable.h"
#include "Foundation/vec3.h"

// -- system includes --
#include <string>
#include <utility>

using std::string;

class MpiWTimers;
class TML_Comm;

/*!
  \brief Abstract base class for sublattices.
*/
class ASubLattice : public esys::lsm::CheckPointable
{
private:
  std::string m_particleType;

protected:

  // -- neighbortable --

public:
  typedef std::pair<int,int>          ParticleIdPair;
  typedef std::vector<ParticleIdPair> ParticleIdPairVector;
  typedef std::vector<int> IdVector;

  virtual ~ASubLattice();
  void setNTSize(int);
  virtual void setParticleType(const std::string &particleType)
  {
    m_particleType = particleType;
  }
  virtual const std::string &getParticleType() const
  {
    return m_particleType;
  }
  /****fluid contents: begin***/
  virtual void addFluidInteraction()=0;
  virtual void updateFluid()=0;  
  virtual void exchangeCells()=0;
  virtual void sendCoeffi()=0;
  virtual void recvPressure()=0;
  virtual void solveMatrix()=0;

  /****fluid contents: end***/

  virtual void setTimeStepSize(double dt) = 0;
  virtual vector<int> getCommCoords() const=0;
  virtual vector<int> getCommDims() const=0;
  virtual void receiveParticles()=0;
  virtual void receiveConnections()=0;
  virtual void addWall()=0;
  virtual void addElasticWIG()=0;
  virtual void addBondedWIG()=0;
  virtual void addDirBondedWIG()=0;
  virtual void addViscWIG()=0;
  virtual void addTaggedElasticWIG()=0;
  virtual void initNeighborTable(const Vec3&,const Vec3&)=0;
  virtual void initNeighborTable(const Vec3&,const Vec3&,const vector<bool>&)=0;
  virtual void addPairIG()=0;
  virtual void addTaggedPairIG()=0;
  virtual void addTriMesh()=0;
  virtual void addTriMeshIG()=0;
  virtual void addBondedTriMeshIG()=0;
  virtual void addMesh2D()=0;
  virtual void addMesh2DIG()=0;
  virtual void addBondedMesh2DIG()=0;
  virtual void addSingleIG()=0;
  virtual void addBondedIG()=0;
  virtual void addCappedBondedIG()=0;
  virtual void addShortBondedIG()=0;
  virtual void addRotBondedIG()=0;
  virtual void addBrittleBeamSCIG()=0;
  virtual void addBrittleBeamDZCIG()=0;
  virtual void addRotThermBondedIG()=0;
  virtual void addDamping()=0;
  //virtual void addRotDamping()=0;
  virtual void setExIG()=0;
  virtual void initComplex();
  virtual void removeIG()=0;
  virtual void getWallPos()=0;
  virtual void getWallForce()=0;
  virtual void addSphereBody()=0;
  virtual void addESphereBodyIG()=0;
  virtual void getSphereBodyPos()=0;
  virtual void getSphereBodyForce()=0;

  virtual const MPI_Comm &getWorkerComm() const = 0;

  virtual void rebuildParticleArray()=0;
  virtual void rebuildInteractions()=0;
  virtual void searchNeighbors()=0;
  virtual void checkNeighbors()=0;

  virtual void updateInteractions()=0;

  virtual int getNumParticles() = 0;

  virtual std::pair<double, int> findParticleNearestTo(const Vec3 &pt) = 0;

  virtual std::pair<int, Vec3> getParticlePosn(int particleId) = 0;

//  virtual ParticleIdPairVector getBondGroupIdPairs(const std::string &groupName) = 0;

  virtual void oneStep()=0;
  virtual void exchangePos()=0;

  // moving stuff around
  virtual void moveParticleTo()=0;
  virtual void changeRadiusBy()=0;
  virtual void moveTaggedParticlesBy() = 0;
  virtual void moveSingleParticleTo(int particleId, const Vec3 &posn)=0;
  virtual void moveWallBy()=0;
  virtual void moveSphereBodyBy()=0;
  virtual void setWallNormal()=0;
  virtual void applyForceToWall()=0;
  virtual void setVelocityOfWall()=0;
  virtual void setParticleVelocity()=0;
  virtual void setParticleDensity()=0;
  virtual void resetParticleRotation()=0;
  virtual void setTaggedParticleVel()=0;
  virtual void setParticleAngularVelocity(){};
  virtual void setParticleNonDynamic()=0;
  virtual void setParticleNonRot()=0;
  virtual void tagParticleNearestTo()=0;
  virtual void moveSingleNode()=0;
  virtual void moveTaggedNodes()=0;
  virtual void translateMeshBy(const std::string &meshName, const Vec3 &translation)=0;
  virtual void rotateMeshBy(const std::string &meshName, const Vec3& origin, const Vec3 &axis, double angle)=0;

  virtual void setTimer(MpiWTimers &timers) = 0;

  virtual void do2dCalculations(bool do2d) = 0;

  // --- setting interaction parameters during a simulation ---
  virtual void setInteractionParameter() = 0;

  // field functions
  virtual void countParticles()=0;

  // "new" field functions
  virtual void addScalarParticleField()=0;
  virtual void addVectorParticleField()=0;
  virtual void addScalarInteractionField()=0;
  virtual void addScalarHistoryInteractionField()=0;
  virtual void addVectorInteractionField()=0;
  virtual void addVectorTriangleField()=0;
  virtual void addScalarTriangleField()=0;
  virtual void sendFieldData()=0;
  virtual void addVectorWallField()=0;
  /****fluid contents: begin***/
  virtual void addScalarFluidField()=0;
  virtual void addVectorFluidField()=0;
  virtual void addScalarFluidInteractionField()=0;
  virtual void addVectorFluidInteractionField()=0;
  /****fluid contents: end***/

  // output
  virtual void printStruct()=0;
  virtual void printData()=0;
  virtual void printTimes()=0;

  // -- mesh data exchange --
  virtual void getMeshNodeRef()=0;
  virtual void getMeshFaceRef()=0;
  virtual void getMesh2DStress()=0;
  virtual void getTriMeshForce()=0;
  virtual void getParticleData(const IdVector &particleIdVector)=0;

  // checkpointing
  virtual void loadCheckPointData(std::istream&){};
};

#endif //__ASUBLATTICE_H
