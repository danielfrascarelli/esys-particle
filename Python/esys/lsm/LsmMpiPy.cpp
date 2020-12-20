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

#include <mpi.h>
#include <boost/version.hpp>
#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/mpl/vector.hpp>
#include "Python/esys/lsm/ParticlePy.h"
#include "Python/esys/lsm/RotParticlePy.h"
#include "Python/esys/lsm/RotParticleViPy.h"
#include "Python/esys/lsm/RotThermalParticlePy.h"
#include "Python/esys/lsm/LsmMpiPy.h"
#include "Parallel/LatticeMaster.h"
#include "Foundation/StringUtil.h"
#include "Python/BoostPythonUtil/ListConverter.h"
#include "Python/BoostPythonUtil/PythonIterIterator.h"
#include "Python/esys/lsm/util/Vec3Py.h"
#include "Python/esys/lsm/util/BoundingBoxPy.h"
#include "Python/esys/lsm/RunnablePy.h"
#include "Python/esys/lsm/geometry/SimpleSpherePy.h"
#include "Python/esys/lsm/geometry/TaggedIdConnectionPy.h"
#include "Python/esys/lsm/CheckPointParamsPy.h"
#include "Python/esys/lsm/InteractionParamsPy.h"
#include "Python/esys/lsm/MeshBuildParamsPy.h"
#include "Python/esys/lsm/BondedTriMeshPrmsPy.h"
#include "Python/esys/lsm/ElasticMesh2DPrmsPy.h"
#include "Python/esys/lsm/ElasticTriMeshPrmsPy.h"
#include "Python/esys/lsm/BondedMesh2DPrmsPy.h"
#include "Python/esys/lsm/WallPrmsPy.h"
#include "Python/esys/lsm/SphereBodyPrmsPy.h"
#include "Python/esys/lsm/ParticleFieldSaverPrmsPy.h"
#include "Python/esys/lsm/InteractionFieldSaverPrmsPy.h"
#include "Python/esys/lsm/WallFieldSaverPrmsPy.h"
#include "Python/esys/lsm/TriangleFieldSaverPrmsPy.h"
#include "Python/esys/lsm/LmParticleAdder.h"
#include "Python/esys/lsm/TriggerPrmsPy.h"
#include "Foundation/console.h"
#include "Python/esys/lsm/FluidFieldSaverPrmsPy.h" //fluid contents
#include "Python/esys/lsm/FluidInteractionFieldSaverPrmsPy.h" //fluid contents

using namespace boost;
using namespace esys::lsm;

#include <vector>
#include <fstream>
#include <stdexcept>
#include <string>

namespace esys
{
  namespace lsm
  {
    class LsmMpiPy::Impl
    {
    public:
      Impl() : m_latticeMaster(), m_nameTypeMap()
      {
      }

      CLatticeMaster m_latticeMaster;
      LsmMpiPy::InteractionNameTypeMap m_nameTypeMap;
    };

    void throwValueError(const std::string &msg)
    {
      PyErr_SetString(PyExc_ValueError, msg.c_str());
      boost::python::throw_error_already_set();
    }

    std::string joinIntVector(
      const std::vector<int> &v,
      const std::string &delim
    )
    {
      return
        StringUtil::join<
          std::vector<int>,
          StringUtil::StdOStreamOp<std::vector<int>::const_iterator>
        >(v, delim);
    }

    void checkParticleTypePy(const std::string &particleType)
    {
      if (
        (particleType != "NRotSphere")
        &&
        (particleType != "RotSphere")
        &&
        (particleType != "RotSphereVi")
        &&
        (particleType != "RotThermalSphere")
      )
      {
        throwValueError(
          std::string()
          +
          "Invalid particle type, needs to be one of"
          +
          " 'NRotSphere', 'RotSphere', 'RotSphereVi' or 'RotThermalSphere'"
          +
          ", got particleType='"
          +
          particleType + "'"
        );
      }
    }

    void checkMpiDimensions(
      int numProcesses,
      const std::vector<int> &mpiDimVector
    )
    {
      if (numProcesses > 0)
      {
        if (mpiDimVector.size() == 3)
        {
          int prod = 1;
          for (unsigned int i = 0; i < mpiDimVector.size(); i++)
          {
            const int dim = mpiDimVector[i];
            if (dim > 0)
            {
              prod *= dim;
            }
            else if (dim < 0)
            {
              throwValueError(
                std::string()
                +
                "MPI dimension list contains negative element, values"
                +
                " must"
                +
                " be >= 0, got MPI dimension list=" + joinIntVector(mpiDimVector, ",")
              );
            }
          }
          if ((numProcesses % prod) != 0)
          {
            throwValueError(
              std::string()
              +
              "Inconsistency beween numProcesses and MPI dimension list,"
              +
              " numProcesses must be a multiple of the product of"
              +
              " non-zero"
              +
              " MPI dimension list elements, got (numProcesses % product)=("
              +
              StringUtil::toString(numProcesses) + " % " + StringUtil::toString(prod)
              +
              ") = " + StringUtil::toString(numProcesses % prod)
            );
          }
        }
        else
        {
          throwValueError(
            "mpiDimList must have 3 elements, got mpiDimList="
            +
            joinIntVector(mpiDimVector, ",")
          );
        }
      }
      else
      {
        throwValueError(
          std::string()
          +
          "Number of processes must be greater than 0, got numProcesses="
          +
          StringUtil::toString(numProcesses)
        );
      }
    }

    void checkMpiDimensionsPy(
      int numWorkerProcesses,
      const boost::python::list &mpiDimList
    )
    {
      checkMpiDimensions(numWorkerProcesses, bpu::listToVector<int>(mpiDimList));
    }

    LsmMpiPy::LsmMpiPy(
      int numWorkerProcesses,
      const python::list &mpiDimList
    )
      : m_implPtr(new Impl)
    {
      getLatticeMaster().setupWorkers(numWorkerProcesses);
      getLatticeMaster().setProcessDims(
        bpu::listToVector<unsigned int>(mpiDimList)
      );
    }

    LsmMpiPy::~LsmMpiPy()
    {
    }

    LsmMpiPy::InteractionNameTypeMap &LsmMpiPy::getNameTypeMap()
    {
      return m_implPtr->m_nameTypeMap;
    }

    const LsmMpiPy::InteractionNameTypeMap &LsmMpiPy::getNameTypeMap() const
    {
      return m_implPtr->m_nameTypeMap;
    }

    int LsmMpiPy::getNumWorkerProcesses() const
    {
      return getLatticeMaster().getNumWorkerProcesses();
    }

    void LsmMpiPy::initVerletModel(
      const std::string &particleType,
      double gridSpacing,
      double verletDist
    )
    {
      checkParticleTypePy(particleType);
      std::string lmParticleType;
      if (particleType == "RotSphere")
      {
        lmParticleType = "Rot";
      }
      else if (particleType == "RotSphereVi")
      {
        lmParticleType = "RotVi";
      }
      else if (particleType == "RotThermalSphere")
      {
        lmParticleType = "RotTherm";
      }
      else
      {
        lmParticleType = "Basic";
      }
      getLatticeMaster().makeLattice(
        lmParticleType.c_str(),
        gridSpacing,
        verletDist
      );
    }


    double LsmMpiPy::getTimeStepSize() const
    {
      return getLatticeMaster().getTimeStepSize();
    }

    void LsmMpiPy::setTimeStepSize(double dt)
    {
      getLatticeMaster().setTimeStepSize(dt);
    }

    std::string LsmMpiPy::getLsmVersion() const
    {
      return getLatticeMaster().getLsmVersion();
    }

    void LsmMpiPy::setTimingFileName(const std::string &fileName)
    {
      getLatticeMaster().setTimingFileName(fileName);
    }

    void LsmMpiPy::setSlaveTimingFileName(const std::string &fileName)
    {
      getLatticeMaster().saveTimingDataToFile(fileName);
    }

    void LsmMpiPy::createConnections(boost::python::object &iteratable)
    {
      bpu::PythonIterIterator<TaggedIdConnectionPy &> it(iteratable);
      getLatticeMaster().addConnections(it);
    }


    /****fluid contents: begin****/
    void LsmMpiPy::createFluidInteraction(
      double cellside,
      double Bw,
      double Bp,
      double Mu,
      double alpha,
      double flowrate,
      double pressure,
      const Vec3Py& inflow,
      const Vec3Py& outflow,
      double dt_f
    )
    {
      getLatticeMaster().addFluidInteraction(cellside,Bw,Bp,Mu,alpha,flowrate,pressure,inflow,outflow,dt_f);
    }

    void LsmMpiPy::createFluidInteractionVec3(
      const Vec3Py& cellside,
      double Bw,
      double Bp,
      double Mu,
      double alpha,
      double flowrate,
      double pressure,
      const Vec3Py& inflow,
      const Vec3Py& outflow,
      double dt_f
    )
    {
      getLatticeMaster().addFluidInteractionVec3(cellside,Bw,Bp,Mu,alpha,flowrate,pressure,inflow,outflow,dt_f);
    }

    void LsmMpiPy::createFluidScalarFieldSaver(
      const FluidScalarFieldSaverPrmsPy &prms
    )
    {
      getLatticeMaster().addScalarFluidSaveField(
        prms.getFileName(),
        prms.getFieldName(),
        prms.getFileFormat(),
        prms.getBeginTimeStep(),
        prms.getEndTimeStep(),
        prms.getTimeStepIncr()
      );
    }

    void LsmMpiPy::createFluidVectorFieldSaver(
      const FluidVectorFieldSaverPrmsPy &prms
    )
    {
      getLatticeMaster().addVectorFluidSaveField(
        prms.getFileName(),
        prms.getFieldName(),
        prms.getFileFormat(),
        prms.getBeginTimeStep(),
        prms.getEndTimeStep(),
        prms.getTimeStepIncr()
      );
    }

    void LsmMpiPy::createFluidInteractionScalarFieldSaver(
      const FluidInteractionScalarFieldSaverPrmsPy &prms
    )
    {
      getLatticeMaster().addScalarFluidInteractionSaveField(
        prms.getFileName(),
        prms.getFieldName(),
        prms.getFileFormat(),
        prms.getBeginTimeStep(),
        prms.getEndTimeStep(),
        prms.getTimeStepIncr()
      );
    }

    void LsmMpiPy::createFluidInteractionVectorFieldSaver(
      const FluidInteractionVectorFieldSaverPrmsPy &prms
    )
    {
      getLatticeMaster().addVectorFluidInteractionSaveField(
        prms.getFileName(),
        prms.getFieldName(),
        prms.getFileFormat(),
        prms.getBeginTimeStep(),
        prms.getEndTimeStep(),
        prms.getTimeStepIncr()
      );
    }
    /****fluid contents: end****/



    void LsmMpiPy::createParticles(boost::python::object &iteratable)
    {
      if (getParticleType() == "NRotSphere")
      {
        LmParticleAdder<
          boost::mpl::vector<RotParticlePy,ParticlePy,SimpleSpherePy>,
          CParticle
        > adder;
        adder.addParticles(iteratable, getLatticeMaster());
      }
      else if (getParticleType() == "RotSphere")
      {
        LmParticleAdder<
          boost::mpl::vector<RotParticlePy,ParticlePy,SimpleSpherePy>,
          CRotParticle
        > adder;
        adder.addParticles(iteratable, getLatticeMaster());
      }
      else if (getParticleType() == "RotSphereVi")
      {
        LmParticleAdder<
          boost::mpl::vector<RotParticleViPy,ParticlePy,SimpleSpherePy>,
          CRotParticleVi
        > adder;
        adder.addParticles(iteratable, getLatticeMaster());
      }
      else if (getParticleType() == "RotThermalSphere")
      {
        LmParticleAdder<
          boost::mpl::vector<RotThermalParticlePy,RotParticleViPy,ParticlePy,SimpleSpherePy>,
          CRotThermParticle
        > adder;
        adder.addParticles(iteratable, getLatticeMaster());
      }
      else
      {
        throw
          std::runtime_error(
            std::string("Unknown particle type:")
            +
            getLatticeMaster().getParticleType()
          );
      }
    }

    void LsmMpiPy::createParticle(boost::python::object &particle)
    {
      boost::python::list l;
      l.append(particle);
      createParticles(l);
    }

    std::string LsmMpiPy::getParticleType() const
    {
      const std::string lmParticleType = getLatticeMaster().getParticleType();
      if (lmParticleType == "Rot")
      {
        return "RotSphere";
      }
      else if (lmParticleType == "RotTherm")
      {
        return "RotThermalSphere";
      }
      else if (lmParticleType == "RotVi")
      {
        return "RotSphereVi";
      }
      else
      {
        return "NRotSphere";
      }
    }

    void LsmMpiPy::readGeometry(const std::string &fileName)
    {
      getLatticeMaster().readGeometryFile(fileName);
    }

    int LsmMpiPy::getNumParticles()
    {
      return getLatticeMaster().getNumParticles();
    }

    int LsmMpiPy::getTimeStep() const
    {
      return getLatticeMaster().getTimeStep();
    }

    // ---  interaction creation functions ---

    void LsmMpiPy::createNRotElasticInteractGrp(
      const NRotElasticPrmsPy &prms
    )
    {
      getNameTypeMap()[prms.getName()] = prms.getTypeString();
      getLatticeMaster().addPairIG(prms);
    }

    void LsmMpiPy::createHertzianElasticIG(
      const HertzianElasticPrmsPy &prms
    )
    {
      getNameTypeMap()[prms.getName()] = prms.getTypeString();
      getLatticeMaster().addPairIG(prms);
    }

    void LsmMpiPy::createHertzianViscoElasticFrictionIG(
      const HertzianViscoElasticFrictionPrmsPy &prms
    )
    {
      HertzianViscoElasticFrictionPrmsPy p = prms;

      p.setTimeStepSize(getTimeStepSize());
      getNameTypeMap()[p.getName()] = p.getTypeString();
      getLatticeMaster().addPairIG(p);
    }

    void LsmMpiPy::createHertzianViscoElasticIG(
      const HertzianViscoElasticPrmsPy &prms
    )
    {
      getNameTypeMap()[prms.getName()] = prms.getTypeString();
      getLatticeMaster().addPairIG(prms);
    }

    void LsmMpiPy::createHertzMindlinIG(
      const HertzMindlinPrmsPy &prms
    )
    {
      HertzMindlinPrmsPy p = prms;

      p.setTimeStepSize(getTimeStepSize());
      getNameTypeMap()[p.getName()] = p.getTypeString();
      getLatticeMaster().addPairIG(p);
    }

    void LsmMpiPy::createHertzMindlinViscoIG(
      const HertzMindlinViscoPrmsPy &prms
    )
    {
      HertzMindlinViscoPrmsPy p = prms;

      p.setTimeStepSize(getTimeStepSize());
      getNameTypeMap()[p.getName()] = p.getTypeString();
      getLatticeMaster().addPairIG(p);
    }

    void LsmMpiPy::createLinearDashpotIG(
      const LinearDashpotPrmsPy &prms
    )
    {
      getNameTypeMap()[prms.getName()] = prms.getTypeString();
      getLatticeMaster().addPairIG(prms);
    }

    void LsmMpiPy::createNRotBondInteractGrp(
      const NRotBondPrmsPy &bondPrms
    )
    {
      getNameTypeMap()[bondPrms.getName()] = bondPrms.getTypeString();
      getLatticeMaster().addBondedIG(bondPrms);
    }

   void LsmMpiPy::createCappedNRotBondInteractGrp(
      const CappedNRotBondPrmsPy &bondPrms
    )
    {
      getNameTypeMap()[bondPrms.getName()] = bondPrms.getTypeString();
      getLatticeMaster().addCappedBondedIG(
        bondPrms.tag,
        bondPrms.getName(),
        bondPrms.k,
        bondPrms.rbreak,
        bondPrms.m_force_limit
      );
    }

    void LsmMpiPy::createNRotShortBondInteractGrp(
      const NRotShortBondPrmsPy &bondPrms
    )
    {
      getNameTypeMap()[bondPrms.getName()] = bondPrms.getTypeString();
      getLatticeMaster().addShortBondedIG(
        bondPrms.tag,
        bondPrms.getName(),
        bondPrms.k,
        bondPrms.rbreak
      );
    }

    void LsmMpiPy::createNRotFrictionInteractGrp(
      const NRotFrictionPrmsPy &prms
    )
    {
      NRotFrictionPrmsPy p = prms;

      p.setTimeStepSize(getTimeStepSize());
      getNameTypeMap()[p.getName()] = p.getTypeString();
      getLatticeMaster().addPairIG(p);
    }

    void LsmMpiPy::createSpringDashpotFrictionInteractGrp(
      const SpringDashpotFrictionPrmsPy &prms
    )
    {
      SpringDashpotFrictionPrmsPy p = prms;

      p.setTimeStepSize(getTimeStepSize());
      getNameTypeMap()[p.getName()] = p.getTypeString();
      getLatticeMaster().addPairIG(p);
    }

    void LsmMpiPy::createRotBondInteractGrp(const RotBondPrmsPy &bondPrms)
    {
      getNameTypeMap()[bondPrms.getName()] = bondPrms.getTypeString();
      getLatticeMaster().addRotBondedIG(
        bondPrms.tag,
        bondPrms.getName(),
        bondPrms.kr,
        bondPrms.ks,
        bondPrms.kt,
        bondPrms.kb,
        bondPrms.max_nForce,
        bondPrms.max_shForce,
        bondPrms.max_tMoment,
        bondPrms.max_bMoment,
        bondPrms.scaling,
        bondPrms.meanR_scaling,
        bondPrms.truncated
      );
    }

    BondInteractionGroupPy LsmMpiPy::createRotThermBondInteractGrp(
      const RotThermBondPrmsPy &bondPrms
    )
    {
      getNameTypeMap()[bondPrms.getName()] = bondPrms.getTypeString();
      getLatticeMaster().addRotThermBondedIG(bondPrms);

      return BondInteractionGroupPy(*this, bondPrms.getName());
    }

    void LsmMpiPy::createBrittleBeamInteractGrp(const BrittleBeamPrmsPy &bondPrms)
    {
      getNameTypeMap()[bondPrms.getName()] = bondPrms.getTypeString();
      getLatticeMaster().addRotBondedIG(
        bondPrms.tag,
        bondPrms.getName(),
        bondPrms.kr,
        bondPrms.ks,
        bondPrms.kt,
        bondPrms.kb,
        bondPrms.max_nForce,
        bondPrms.max_shForce,
        bondPrms.max_tMoment,
        bondPrms.max_bMoment,
        bondPrms.scaling,
        bondPrms.meanR_scaling,
        bondPrms.truncated
      );
    }

    void LsmMpiPy::createBrittleBeamSCInteractGrp(const BrittleBeamSCPrmsPy &bondPrms)
    {
        getNameTypeMap()[bondPrms.getName()] = bondPrms.getTypeString();
        getLatticeMaster().addBrittleBeamSCIG(bondPrms);
    }

    void LsmMpiPy::createBrittleBeamDZCInteractGrp(const BrittleBeamDZCPrmsPy &bondPrms)
    {
        getNameTypeMap()[bondPrms.getName()] = bondPrms.getTypeString();
        getLatticeMaster().addBrittleBeamDZCIG(bondPrms);
    }
    
    void LsmMpiPy::createFrictionInteractGrp(const FrictionPrmsPy &prms)
    {
      FrictionPrmsPy p=prms;

      p.setTimeStepSize(getTimeStepSize());
      getNameTypeMap()[p.getName()] = p.getTypeString();
      getLatticeMaster().addPairIG(p);
    }

    void LsmMpiPy::createRotFrictionInteractGrp(const RotFrictionPrmsPy &prms)
    {
      RotFrictionPrmsPy p=prms;

      p.setTimeStepSize(getTimeStepSize());
      getNameTypeMap()[p.getName()] = p.getTypeString();
      getLatticeMaster().addPairIG(p);
    }

    void LsmMpiPy::createRotThermFrictionInteractGrp(
      const RotThermFrictionPrmsPy &prms
    )
    {
      getNameTypeMap()[prms.getName()] = prms.getTypeString();
      getLatticeMaster().addPairIG(prms);
    }

    void LsmMpiPy::createVWFrictionIG(const VWFrictionPrmsPy &prms)
    {
      VWFrictionPrmsPy p=prms;

      p.setTimeStepSize(getTimeStepSize());
      getNameTypeMap()[prms.getName()] = prms.getTypeString();
      getLatticeMaster().addPairIG(prms);
    }

    void LsmMpiPy::createRotElasticInteractGrp(const RotElasticPrmsPy& prms)
    {
      getNameTypeMap()[prms.getName()] = prms.getTypeString();
      getLatticeMaster().addPairIG(prms);
    }

    void LsmMpiPy::createRotThermElasticInteractGrp(
      const RotThermElasticPrmsPy &prms
    )
    {
      getNameTypeMap()[prms.getName()] = prms.getTypeString();
      getLatticeMaster().addPairIG(prms);
    }

    void LsmMpiPy::createDamping(const DampingPrmsPy &prms)
    {
      DampingPrmsPy p = prms;
      p.setTimeStepSize(getTimeStepSize());
      getLatticeMaster().addDamping(p);
    }

    void LsmMpiPy::createLocalDamping(const LocalDampingPrmsPy &prms)
    {
      LocalDampingPrmsPy p = prms;
      p.setTimeStepSize(getTimeStepSize());
      getLatticeMaster().addDamping(p);
    }

    void LsmMpiPy::createRotLocalDamping(const RotLocalDampingPrmsPy &prms)
    {
      RotLocalDampingPrmsPy p = prms;
      p.setTimeStepSize(getTimeStepSize());
      getLatticeMaster().addDamping(p);
    }

    void LsmMpiPy::createABCDamping(const ABCDampingPrmsPy &prms)
    {
      ABCDampingPrmsPy p = prms;
      p.setTimeStepSize(getTimeStepSize());
      getLatticeMaster().addDamping(p);
    }

    void LsmMpiPy::createGravity(const GravityPrmsPy& prms)
    {
      getLatticeMaster().addSingleIG(prms);
    }

    void LsmMpiPy::createBuoyancy(const BuoyancyPrmsPy& prms)
    {
      getLatticeMaster().addSingleIG(prms);
    }

    void  LsmMpiPy::removeInteractionGrp(const std::string& name)
    {
      getLatticeMaster().removeIG(name);
    }

    // ----------------------------------------------
    //     tagged  interaction creation functions
    // ----------------------------------------------

    void LsmMpiPy::createRotFrictionInteractGrpTag(const RotFrictionPrmsPy &prms,
                           int tag1, int mask1,
                           int tag2, int mask2)
    {
      RotFrictionPrmsPy p=prms;

      p.setTimeStepSize(getTimeStepSize());
      getNameTypeMap()[p.getName()] = p.getTypeString();
      getLatticeMaster().addTaggedPairIG(p,tag1,mask1,tag2,mask2);
    }

    void LsmMpiPy::createFrictionInteractGrpTag(const FrictionPrmsPy &prms,
                                                   int tag1, int mask1,
                                                   int tag2, int mask2)
    {
      FrictionPrmsPy p=prms;

      p.setTimeStepSize(getTimeStepSize());
      getNameTypeMap()[p.getName()] = p.getTypeString();
      getLatticeMaster().addTaggedPairIG(p,tag1,mask1,tag2,mask2);
    }


    void LsmMpiPy::createNRotFrictionInteractGrpTag(const NRotFrictionPrmsPy &prms,
                           int tag1, int mask1,
                           int tag2, int mask2)
    {
      NRotFrictionPrmsPy p=prms;

      p.setTimeStepSize(getTimeStepSize());
      getNameTypeMap()[p.getName()] = p.getTypeString();
      std::cerr << "createNRotFrictionInteractGrpTag " << tag1 << " , "<< mask1 << " , "<< tag2 << " , "<< mask2 << " , " << std::endl;
      getLatticeMaster().addTaggedPairIG(p,tag1,mask1,tag2,mask2);
    }

    void LsmMpiPy::createSpringDashpotFrictionInteractGrpTag(const SpringDashpotFrictionPrmsPy &prms,
                                                   int tag1, int mask1,
                                                   int tag2, int mask2)
    {
      SpringDashpotFrictionPrmsPy p=prms;

      p.setTimeStepSize(getTimeStepSize());
      getNameTypeMap()[p.getName()] = p.getTypeString();
      std::cerr << "createNRotFrictionInteractGrpTag " << tag1 << " , "<< mask1 << " , "<< tag2 << " , "<< mask2 << " , " << std::endl;
      getLatticeMaster().addTaggedPairIG(p,tag1,mask1,tag2,mask2);
    }

    void LsmMpiPy::createLinearDashpotInteractGrpTag(const LinearDashpotPrmsPy &prms,
                             int tag1, int mask1,
                             int tag2, int mask2)
    {
      LinearDashpotPrmsPy p=prms;

      getNameTypeMap()[p.getName()] = p.getTypeString();
      getLatticeMaster().addTaggedPairIG(p,tag1,mask1,tag2,mask2);
    }

    void LsmMpiPy::createElasticInteractGrpTag(const NRotElasticPrmsPy &prms,
      int tag1, int mask1,
      int tag2, int mask2)
    {
      NRotElasticPrmsPy p=prms;

      getNameTypeMap()[p.getName()] = p.getTypeString();
      getLatticeMaster().addTaggedPairIG(p,tag1,mask1,tag2,mask2);
    }

    void LsmMpiPy::createRotElasticInteractGrpTag(
      const RotElasticPrmsPy& prms,
      int tag1,
      int mask1,
      int tag2,
      int mask2
    ) {
      RotElasticPrmsPy p=prms;

      getNameTypeMap()[p.getName()] = p.getTypeString();
      getLatticeMaster().addTaggedPairIG(p,tag1,mask1,tag2,mask2);
    }

    // ----------------------------------------------
    //     interaction group exclusion
    // ----------------------------------------------

    void LsmMpiPy::createExclusion(
      const std::string &interactionName1,
      const std::string &interactionName2
    )
    {
      getLatticeMaster().addExIG(interactionName1, interactionName2);
    }

    // -- checkpoint & snapshot setup ---
    // restart checkpoints
    void LsmMpiPy::createCheckPointer(const RestartCheckPointPrmsPy &prms)
    {
      getLatticeMaster().performCheckPoints(
        prms.getFileNamePrefix(),
        prms.getBeginTimeStep(),
        prms.getEndTimeStep(),
        prms.getTimeStepIncr(),
	prms.getPrecision()
      );
    }

    void LsmMpiPy::createCheckPointerThroughMaster(const RestartCheckPointPrmsPy &prms)
    {
      getLatticeMaster().performCheckPointsThroughMaster(
        prms.getFileNamePrefix(),
        prms.getBeginTimeStep(),
        prms.getEndTimeStep(),
        prms.getTimeStepIncr(),
	prms.getPrecision()
      );
    }

    // visualisation snapshots
    void LsmMpiPy::createSnapShots(const CheckPointPrmsPy &prms)
    {
      getLatticeMaster().initSnapShotController(
        prms.getFileNamePrefix(),
        prms.getBeginTimeStep(),
        prms.getEndTimeStep(),
        prms.getTimeStepIncr()
      );
    }

    void LsmMpiPy::loadCheckPoint(const std::string& filename)
    {
      getLatticeMaster().loadCheckPointData(filename);
    }

    int LsmMpiPy::getNumTimeSteps() const
    {
      return getLatticeMaster().getNumSteps();
    }

    void LsmMpiPy::setNumTimeSteps(int numTimeSteps)
    {
      getLatticeMaster().setNumSteps(numTimeSteps);
    }

    //-------------------
    // Mesh related functions
    //-------------------
    /*!
      read triangle mesh from file
    */
    void LsmMpiPy::readMeshWithTag(const std::string &fileName, const std::string &meshName,int tag)
    {
      getLatticeMaster().readAndDistributeTriMesh(meshName, fileName, tag);
    }

    /*!
      read triangle mesh from file ignoring tag
    */
    void LsmMpiPy::readMesh(const std::string &fileName, const std::string &meshName)
    {
      getLatticeMaster().readAndDistributeTriMesh(meshName, fileName);
    }

    /*
       create triangle mesh from within the python script
    */
    void LsmMpiPy::createTriMesh(
      const std::string &meshName,
      const boost::python::object &nodeSequence,
      const boost::python::object &triSequence
    )
    {
      const std::size_t numNodes = bpu::len(nodeSequence);
      MeshNodeDataVector nodeVector;
      for (std::size_t i = 0; i < numNodes; i++)
      {
        boost::python::object node = nodeSequence[i];
        int tag = 0;
        if (bpu::len(node) > 2)
        {
          tag = boost::python::extract<int>(node[2]);
        }
        nodeVector.push_back(
          MeshNodeData(
            boost::python::extract<int>(node[0]),
            Vec3Py(node[1]),
            tag
          )
        );
      }
      const std::size_t numTriangles = bpu::len(triSequence);
      MeshTriDataVector triVector;
      for (std::size_t i = 0; i < numTriangles; i++)
      {
        boost::python::object triangle = triSequence[i];
        int tag = 0;
        if (bpu::len(triangle) > 2)
        {
          tag = boost::python::extract<int>(triangle[2]);
        }
        triVector.push_back(
          MeshTriData(
            boost::python::extract<int>(triangle[0]),
            boost::python::extract<int>(triangle[1][0]),
            boost::python::extract<int>(triangle[1][1]),
            boost::python::extract<int>(triangle[1][2]),
            tag
          )
        );
      }
      getLatticeMaster().createTriMesh(meshName, nodeVector, triVector);
    }

    void LsmMpiPy::readMesh2D(
      const std::string &fileName,
      const std::string &meshName,
      int tag
    )
    {
      getLatticeMaster().addMesh2D(meshName, fileName,tag);
    }

    void LsmMpiPy::translateMesh(
      const std::string &meshName,
      const Vec3Py &translation
    )
    {
      getLatticeMaster().translateMeshBy(meshName, translation);
    }

    void LsmMpiPy::rotateMesh(
      const std::string &meshName,
      const Vec3Py& origin,
      const Vec3Py& axis,
      const double angle
    )
    {
      getLatticeMaster().rotateMeshBy(meshName, origin, axis, angle);
    }

    void LsmMpiPy::createNRotElasticTriMeshInteractGrp(
      const NRotElasticTriMeshPrmsPy &prms
    )
    {
      getLatticeMaster().addTriMeshIG(prms);
    }

    void LsmMpiPy::createNRotBondedTriMeshInteractGrp(
      const NRotBondedTriMeshPrmsPy &prms
    )
    {
      if (prms.haveTagBuildPrms()) {
        getLatticeMaster().addBondedTriMeshIG(prms, prms.getTagBuildPrms());
      }
      else if (prms.haveGapBuildPrms()) {
        getLatticeMaster().addBondedTriMeshIG(prms, prms.getGapBuildPrms());
      }
      else {
        throw std::runtime_error("Unknown bonded triangular mesh build prms.");
      }
    }

    void LsmMpiPy::createNRotElasticMesh2DInteractGrp(
      const NRotElasticMesh2DPrmsPy &prms
    )
    {
      getLatticeMaster().addMesh2DIG(prms);
    }

    void LsmMpiPy::createNRotElasticLinMeshInteractGrp(
      const NRotElasticLinMeshPrmsPy &prms
    )
    {
      getLatticeMaster().addMesh2DIG(prms);
    }

    void LsmMpiPy::createNRotBondedLinMeshInteractGrp(
      const NRotBondedLinMeshPrmsPy &prms
    )
    {
      if (prms.haveTagBuildPrms()) {
        getLatticeMaster().addBondedMesh2DIG(prms, prms.getTagBuildPrms());
      }
      else if (prms.haveGapBuildPrms()) {
        getLatticeMaster().addBondedMesh2DIG(prms, prms.getGapBuildPrms());
      }
      else {
        throw std::runtime_error("Unknown bonded triangular mesh build prms.");
      }
    }

    void LsmMpiPy::moveSingleMeshNodeBy(
      const string& meshname,
      int id,
      const Vec3Py& d
    )
    {
      getLatticeMaster().moveSingleNodeBy(meshname,id,d);
    }

    void LsmMpiPy::addPreTimeStepRunnable(RunnablePy &runnable)
    {
      getLatticeMaster().addPreTimeStepRunnable(runnable);
    }

    void LsmMpiPy::addPostTimeStepRunnable(RunnablePy &runnable)
    {
      getLatticeMaster().addPostTimeStepRunnable(runnable);
    }

    void LsmMpiPy::runTimeStep()
    {
      getLatticeMaster().runOneStep();
    }

    void LsmMpiPy::run()
    {
      getLatticeMaster().run();
    }

     // Exit the simulation after running a series of single steps
     // of the time-integration method.
   void LsmMpiPy::exit()
    {
      getLatticeMaster().runEnd();
    }

    void LsmMpiPy::force2dComputations(bool do2d)
    {
      getLatticeMaster().do2dCalculations(do2d);
    }

    void LsmMpiPy::setBBoxSpatialDomain(const BoundingBoxPy &domain)
    {
      getLatticeMaster().setSpatialDomain(
        domain.getMinPt(),
        domain.getMaxPt()
      );
    }

    void LsmMpiPy::setBBoxSpatialDomainWithCirc(
      const BoundingBoxPy &domain,
      const boost::python::list &circDimList
    )
    {
      getLatticeMaster().setSpatialDomain(
        domain.getMinPt(),
        domain.getMaxPt(),
        bpu::listToVector<int>(circDimList)
      );
    }

    void LsmMpiPy::setSpatialDomain(const Vec3Py &minPt, const Vec3Py &maxPt)
    {
      getLatticeMaster().setSpatialDomain(minPt, maxPt);
    }

    int LsmMpiPy::findClosestParticle(const Vec3Py &pt)
    {
      return getLatticeMaster().findParticleNearestTo(pt);
    }

    Vec3Py LsmMpiPy::getParticlePosn(int particleId)
    {
      return Vec3Py(getLatticeMaster().getParticlePosn(particleId));
    }

    void LsmMpiPy::changeRadiusBy (int tag, const double &deltaR)
    {
       getLatticeMaster().changeRadiusBy(tag, deltaR);
    }

    // -----------------------
    //  moving particles
    // -----------------------
    void LsmMpiPy::moveTaggedParticlesTo(int tag, const Vec3Py &pt)
    {
      getLatticeMaster().moveParticleTo(tag, pt);
    }


    void LsmMpiPy::moveTaggedParticlesBy(int tag, const Vec3Py &displacement)
    {
      getLatticeMaster().moveTaggedParticlesBy(tag, displacement);
    }

    void LsmMpiPy::moveSingleParticleTo(int particleId, const Vec3Py &pt)
    {
      getLatticeMaster().moveSingleParticleTo(particleId, pt);
    }

    // -----------------------
    // Wall related functions
    // -----------------------
    void LsmMpiPy::createWall(
      const string& name,
      const Vec3Py& pos, const Vec3Py& norm
    )
    {
      getLatticeMaster().addWall(name,pos,norm);
    }

    void LsmMpiPy::createSphereBody(
      const string& name,
      const Vec3Py& pos, const double& radius
    )
    {
      getLatticeMaster().addSphereBody(name,pos,radius);
    }

    void LsmMpiPy::moveWallBy(const string& name, const Vec3Py &disp)
    {
      getLatticeMaster().moveWallBy(name,disp);
    }

    void LsmMpiPy::moveSphereBodyBy(const string& name, const Vec3Py &disp)
    {
      getLatticeMaster().moveSphereBodyBy(name,disp);
    }

    void LsmMpiPy::setWallNormal(const string& name, const Vec3Py &wn)
    {
      getLatticeMaster().setWallNormal(name,wn);
    }

    void LsmMpiPy::createNRotBondedWall(
      const NRotBondedWallPrmsPy &prms
    )
    {
      getLatticeMaster().addWallIG(prms);
    }

    void LsmMpiPy::createNRotElasticWall(
      const NRotElasticWallPrmsPy &prms
    )
    {
      getLatticeMaster().addWallIG(prms);
    }

    void LsmMpiPy::createNRotElasticSphereBody(
      const NRotElasticSphereBodyPrmsPy &prms
    )
    {
      getLatticeMaster().addSphereBodyIG(prms);
    }

     void LsmMpiPy::createNRotSoftBondedWall(const NRotSoftBondedWallPrmsPy &prms)
     {
       getLatticeMaster().addWallIG(prms);
     }


    void LsmMpiPy::createNRotElasticWallTagged(const NRotElasticWallPrmsPy &prms, int tag, int mask)
    {
      getLatticeMaster().addTaggedWallIG(prms,tag,mask);
    }

    void LsmMpiPy::applyForceToWall(const string& name, const Vec3Py &Frc)
    {
      getLatticeMaster().applyForceToWall(name,Frc);
    }

    Vec3Py LsmMpiPy::getWallPosition(const std::string& name)
    {
      return Vec3Py(getLatticeMaster().getWallPosn(name));
    }

    Vec3Py LsmMpiPy::getSphereBodyPosition(const std::string& name)
    {
      return Vec3Py(getLatticeMaster().getSphereBodyPosn(name));
    }

    Vec3Py LsmMpiPy::getWallForce(const std::string& name)
    {
      return Vec3Py(getLatticeMaster().getWallForce(name));
    }

    Vec3Py LsmMpiPy::getSphereBodyForce(const std::string& name)
    {
      return Vec3Py(getLatticeMaster().getSphereBodyForce(name));
    }

    // --- particle property setting functions ---
    void LsmMpiPy::setParticleVel(int id,const Vec3Py &V)
    {
      getLatticeMaster().setParticleVel(id,V);
    }

    void LsmMpiPy::setParticleAngVel(int id,const Vec3Py &AV)
    {
      getLatticeMaster().setParticleAngVel(id,AV);
    }

    void LsmMpiPy::setParticleDensity(int tag,int mask, double rho)
    {
      getLatticeMaster().setParticleDensity(tag,mask,rho);
    }

    void LsmMpiPy::resetParticleOrientation(int tag,int mask)
    {
	getLatticeMaster().resetParticleRotation(tag,mask);
    }
    
    void LsmMpiPy::setTaggedParticleVel(int tag, const Vec3Py &V)
    {
      getLatticeMaster().setTaggedParticleVel(tag,V);
    }

    void LsmMpiPy::setVelocityOfWall(const std::string &name,const Vec3Py &V)
    {
      getLatticeMaster().setVelocityOfWall(name,V);
    }

    void LsmMpiPy::tagParticleNearestTo(int tag,int mask,const Vec3Py &Pos)
    {
      getLatticeMaster().tagParticleNearestTo(tag,mask,Pos);
    }

    void LsmMpiPy::setParticleNonDynamic(int tag)
    {
      getLatticeMaster().setParticleNonDynamic(tag);
    }

    void LsmMpiPy::setParticleNonRot(int tag)
    {
      getLatticeMaster().setParticleNonRot(tag);
    }

    void LsmMpiPy::setParticleNonTrans(int id)
    {
      getLatticeMaster().setParticleNonTrans(id);
    }

    // --- modifying interaction parameters ---
    void LsmMpiPy::setInteractionParameter(const std::string& igname,const std::string& pname,double val)
    {
        getLatticeMaster().setInteractionParameter(igname,pname,val);
    }

    // ------------------------------
    //    field saving functions
    // ------------------------------

    void LsmMpiPy::createParticleScalarFieldSaver(
      const ParticleScalarFieldSaverPrmsPy &prms
    )
    {
      getLatticeMaster().addScalarParticleSaveField(
        prms.getFileName(),
        prms.getFieldName(),
        prms.getFileFormat(),
        prms.getBeginTimeStep(),
        prms.getEndTimeStep(),
        prms.getTimeStepIncr()
      );
    }

    void LsmMpiPy::createParticleVectorFieldSaver(const ParticleVectorFieldSaverPrmsPy &prms)
    {
      getLatticeMaster().addVectorParticleSaveField(
        prms.getFileName(),
        prms.getFieldName(),
        prms.getFileFormat(),
        prms.getBeginTimeStep(),
        prms.getEndTimeStep(),
        prms.getTimeStepIncr()
      );
    }

    void LsmMpiPy::createInteractionScalarFieldSaver(
      const InteractionScalarFieldSaverPrmsPy &prms
    )
    {
      if (
        getNameTypeMap().find(prms.getInteractionName())
        !=
        getNameTypeMap().end()
      )
      {
        getLatticeMaster().addScalarInteractionSaveField(
          prms.getFileName(),
          prms.getFieldName(),
          getNameTypeMap()[prms.getInteractionName()],
          prms.getInteractionName(),
          prms.getFileFormat(),
          prms.getBeginTimeStep(),
          prms.getEndTimeStep(),
          prms.getTimeStepIncr()
        );
      }
      else
      {
        std::stringstream msg;
        msg
          << "No interaction named '" << prms.getInteractionName()
          << "' has been created in this model.";
        throw std::runtime_error(msg.str().c_str());
      }
    }

    void LsmMpiPy::createCheckedInteractionScalarFieldSaver(
      const CheckedInteractionScalarFieldSaverPrmsPy &prms
    )
    {
      if (
        getNameTypeMap().find(prms.getInteractionName())
        !=
        getNameTypeMap().end()
      )
      {
        getLatticeMaster().addScalarInteractionSaveField(
          prms.getFileName(),
          prms.getFieldName(),
          getNameTypeMap()[prms.getInteractionName()],
          prms.getInteractionName(),
          prms.getFileFormat(),
          prms.getBeginTimeStep(),
          prms.getEndTimeStep(),
          prms.getTimeStepIncr(),
	  true
        );
      }
      else
      {
        std::stringstream msg;
        msg
          << "No interaction named '" << prms.getInteractionName()
          << "' has been created in this model.";
        throw std::runtime_error(msg.str().c_str());
      }
    }

   void LsmMpiPy::createTaggedInteractionScalarFieldSaver(
      const TaggedInteractionScalarFieldSaverPrmsPy &prms
    )
    {
      if (
        getNameTypeMap().find(prms.getInteractionName())
        !=
        getNameTypeMap().end()
      )
      {
        getLatticeMaster().addTaggedScalarInteractionSaveField(
          prms.getFileName(),
          prms.getFieldName(),
          getNameTypeMap()[prms.getInteractionName()],
          prms.getInteractionName(),
          prms.getFileFormat(),
          prms.getBeginTimeStep(),
          prms.getEndTimeStep(),
          prms.getTimeStepIncr(),
	  prms.getTag(),
	  prms.getMask(),
	  false
        );
      }
      else
      {
        std::stringstream msg;
        msg
          << "No interaction named '" << prms.getInteractionName()
          << "' has been created in this model.";
        throw std::runtime_error(msg.str().c_str());
      }
    }



    void LsmMpiPy::createInteractionVectorFieldSaver(
      const InteractionVectorFieldSaverPrmsPy &prms
    )
    {
      if (
        getNameTypeMap().find(prms.getInteractionName())
        !=
        getNameTypeMap().end()
      )
      {
        getLatticeMaster().addVectorInteractionSaveField(
          prms.getFileName(),
          prms.getFieldName(),
          getNameTypeMap()[prms.getInteractionName()],
          prms.getInteractionName(),
          prms.getFileFormat(),
          prms.getBeginTimeStep(),
          prms.getEndTimeStep(),
          prms.getTimeStepIncr()
        );
      }
      else
      {
        std::stringstream msg;
        msg
          << "No interaction named '" << prms.getInteractionName()
          << "' has been created in this model.";
        throw std::runtime_error(msg.str().c_str());
      }
    }


    void LsmMpiPy::createCheckedInteractionVectorFieldSaver(
      const CheckedInteractionVectorFieldSaverPrmsPy &prms
    )
    {
      if (
        getNameTypeMap().find(prms.getInteractionName())
        !=
        getNameTypeMap().end()
      )
      {
        getLatticeMaster().addVectorInteractionSaveField(
          prms.getFileName(),
          prms.getFieldName(),
          getNameTypeMap()[prms.getInteractionName()],
          prms.getInteractionName(),
          prms.getFileFormat(),
          prms.getBeginTimeStep(),
          prms.getEndTimeStep(),
          prms.getTimeStepIncr(),
	  true
        );
      }
      else
      {
        std::stringstream msg;
        msg
          << "No interaction named '" << prms.getInteractionName()
          << "' has been created in this model.";
        throw std::runtime_error(msg.str().c_str());
      }
    }

    void LsmMpiPy::createTaggedParticleScalarFieldSaver(const TaggedParticleScalarFieldSaverPrmsPy& prms)
    {
      getLatticeMaster().addTaggedScalarParticleSaveField(
	prms.getFileName(),
        prms.getFieldName(),
        prms.getFileFormat(),
        prms.getBeginTimeStep(),
        prms.getEndTimeStep(),
        prms.getTimeStepIncr(),
	prms.getTag(),
	prms.getMask());
    }

    void LsmMpiPy::createTaggedParticleVectorFieldSaver(const TaggedParticleVectorFieldSaverPrmsPy& prms)
    {
      getLatticeMaster().addTaggedVectorParticleSaveField(
	prms.getFileName(),
        prms.getFieldName(),
        prms.getFileFormat(),
        prms.getBeginTimeStep(),
        prms.getEndTimeStep(),
        prms.getTimeStepIncr(),
	prms.getTag(),
	prms.getMask());
    }

    void LsmMpiPy::addTaggedScalarParticleDistributionSaver(
      const string &filename,
      const string& fieldname,
      const string& savetype,
      int t_0,int t_end,int dt,int t_snap,
      int tag,int mask,
      double x_0,double x_max,int nx
    )
    {
      getLatticeMaster().addTaggedScalarParticleDistributionSaver(
        filename,fieldname,
        savetype,t_0,t_end,dt,t_snap,tag,mask,x_0,x_max,nx
      );
    }

    void LsmMpiPy::addVectorTriangleSaveField(const TriangleVectorFieldSaverPrmsPy& prms)
    {
      getLatticeMaster().addVectorTriangleSaveField(
        prms.getFileName(),
        prms.getFieldName(),
	prms.getMeshName(),
        prms.getFileFormat(),
        prms.getBeginTimeStep(),
        prms.getEndTimeStep(),
        prms.getTimeStepIncr()
      );
    }

    void LsmMpiPy::addScalarTriangleSaveField(const TriangleScalarFieldSaverPrmsPy& prms)
    {
      getLatticeMaster().addScalarTriangleSaveField(
        prms.getFileName(),
        prms.getFieldName(),
	prms.getMeshName(),
        prms.getFileFormat(),
        prms.getBeginTimeStep(),
        prms.getEndTimeStep(),
        prms.getTimeStepIncr()
      );
    }

    /*!
      wrap LatticeMaster::addVectorWallField
    */
    void LsmMpiPy::addVectorWallField(
      const WallVectorFieldSaverPrmsPy &prms
    )
    {
      getLatticeMaster().addVectorWallField(
        prms.getFileName(),
        prms.getFieldName(),
        prms.getWallNameVector(),
        prms.getFileFormat(),
        prms.getBeginTimeStep(),
        prms.getEndTimeStep(),
        prms.getTimeStepIncr()
      );
    }

    // --- with trigger ---

    /*!
      wrap
    */
    void  LsmMpiPy::createParticleVectorFieldSaverWithTrigger(
       const MaxTriggerPrmsPy &tprms,
       const ParticleVectorFieldSaverPrmsPy &prms)
    {
      const MaxTrigParams mtp=MaxTrigParams(tprms);

      getLatticeMaster().addVectorParticleSaveFieldWT(
	prms.getFileName(),
        prms.getFieldName(),
        prms.getFileFormat(),
        prms.getBeginTimeStep(),
        prms.getEndTimeStep(),
        prms.getTimeStepIncr(),
	mtp
      );
    }

    /*!
      wrap
    */
    void LsmMpiPy::createTaggedParticleVectorFieldSaverWithTrigger(
       const MaxTriggerPrmsPy &tprms,
       const TaggedParticleVectorFieldSaverPrmsPy& prms)
    {
      const MaxTrigParams mtp=MaxTrigParams(tprms);

      getLatticeMaster().addTaggedVectorParticleSaveFieldWT(
	prms.getFileName(),
        prms.getFieldName(),
        prms.getFileFormat(),
        prms.getBeginTimeStep(),
        prms.getEndTimeStep(),
        prms.getTimeStepIncr(),
	prms.getTag(),
	prms.getMask(),
	mtp
      );
    }

    // --- CONSOLE FUNCTIONS ---

    /*!
      set verbosity via LatticeMaster - call gets communicated to Workers
      boolean argument: false -> 0, true -> 7
    */
    void LsmMpiPy::SetVerbosityPy(bool verbose)
    {
    if(verbose){
        getLatticeMaster().setVerbosity(7);
    } else {
        getLatticeMaster().setVerbosity(0);
    }
    }

    /*!
      set verbosity via LatticeMaster - call gets communicated to Workers
      integer argument
    */
    void LsmMpiPy::SetVerbosityLevelPy(int verbose)
    {
    getLatticeMaster().setVerbosity(verbose);
    }

    /*!
      set "base" filename of the console output files
    */
    void LsmMpiPy::SetConsoleFilenamePy(const std::string& fname)
    {
        getLatticeMaster().setConsoleFilename(fname);
    }


    /*!
      set buffering mode and buffer size of the console
    */
    void LsmMpiPy::SetConsoleBufferedPy(unsigned int bsize)
    {
        getLatticeMaster().setConsoleBuffered(bsize);
    }

    // ------

    const CLatticeMaster &LsmMpiPy::getLatticeMaster() const
    {
      return m_implPtr->m_latticeMaster;
    }

    CLatticeMaster &LsmMpiPy::getLatticeMaster()
    {
      return m_implPtr->m_latticeMaster;
    }

    /*!
      Global setVerbosity call. Only influences the Master
      boolean argument: false -> 0, true -> 7
    */
    void setVerbosityPy(bool verbose)
    {
    if(verbose){
        console.SetVerbose(7);
    } else {
        console.SetVerbose(0);
    }
    }

    /*!
      Global setVerbosity call. Only influences the Master
      integer argument
    */
    void setVerbosityLevelPy(int verbose)
    {
    console.SetVerbose(verbose);
    }

    class NodeRefVisitor
    {
    public:
      NodeRefVisitor(boost::python::object pyObject) : m_pyObject(pyObject)
      {
      }

      void visitNodeRef(int nodeRef)
      {
        m_pyObject.attr("visitNodeRef")(nodeRef);
      }

      private:
        boost::python::object m_pyObject;
    };

    class RefStressVisitor
    {
    public:
      RefStressVisitor(boost::python::object pyObject) : m_pyObject(pyObject)
      {
      }

      void visitRefStressPair(int nodeRef, const Vec3 &force)
      {
        m_pyObject.attr("visitRefStressPair")(nodeRef, Vec3Py(force));
      }

      private:
        boost::python::object m_pyObject;
    };

    void LsmMpiPy::visitNodeRefs2d(const std::string &meshName, boost::python::object pyObject)
    {
      NodeRefVisitor visitor(pyObject);
      getLatticeMaster().visitMesh2dNodeReferences(meshName, visitor);
    }

    void LsmMpiPy::visitRefStressPairs2d(
      const std::string &meshName,
      boost::python::object pyObject
    )
    {
      RefStressVisitor visitor(pyObject);
      getLatticeMaster().visitMesh2dEdgeStress(meshName, visitor);
    }

    void LsmMpiPy::updateInteractions()
    {
      // getLatticeMaster().updateInteractions();
    }

    class ParticleVisitor
    {
    public:
      ParticleVisitor(boost::python::object pyObject) : m_pyObject(pyObject)
      {
      }

      void visitParticle(const CParticle &particle)
      {
        m_pyObject.attr("visitNRotSphere")(ParticlePy(particle));
      }

      void visitRotParticle(const CRotParticle &particle)
      {
        m_pyObject.attr("visitRotSphere")(RotParticlePy(particle));
      }

      void visitRotParticleVi(const CRotParticleVi &particle)
      {
        m_pyObject.attr("visitRotSphereVi")(RotParticleViPy(particle));
      }

      void visitRotThermParticle(const CRotThermParticle &particle)
      {
        m_pyObject.attr("visitRotThermalSphere")(RotThermalParticlePy(particle));
      }

      private:
        boost::python::object m_pyObject;
    };

    class ListGatherVisitor
    {
    public:
      ListGatherVisitor() : m_pyList()
      {
      }

      void visitParticle(const CParticle &particle)
      {
        m_pyList.append(ParticlePy(particle));
      }

      void visitRotParticle(const CRotParticle &particle)
      {
        m_pyList.append(RotParticlePy(particle));
      }

      void visitRotParticleVi(const CRotParticleVi &particle)
      {
        m_pyList.append(RotParticleViPy(particle));
      }

      void visitRotThermParticle(const CRotThermParticle &particle)
      {
        m_pyList.append(RotThermalParticlePy(particle));
      }

      const boost::python::list &getList() const
      {
        return m_pyList;
      }

      private:
        boost::python::list m_pyList;
    };

    class BBoxGatherVisitor
    {
    public:
      BBoxGatherVisitor(const BoundingBoxPy &bbox) : m_pyList(),m_bbox(bbox)
      {
      }

      void visitParticle(const CParticle &particle)
      {
        if (insideBoundingBox(particle)) {
           m_pyList.append(ParticlePy(particle));
        }
      }

      void visitRotParticle(const CRotParticle &particle)
      {
        if (insideBoundingBox(particle)) {
           m_pyList.append(RotParticlePy(particle));
        }
      }

      void visitRotParticleVi(const CRotParticleVi &particle)
      {
        if (insideBoundingBox(particle)) {
           m_pyList.append(RotParticleViPy(particle));
        }
      }

      void visitRotThermParticle(const CRotThermParticle &particle)
      {
        if (insideBoundingBox(particle)) {
           m_pyList.append(RotThermalParticlePy(particle));
        }
      }

      const boost::python::list &getList() const
      {
        return m_pyList;
      }

    private:
      template <typename TmplParticle>
      bool insideBoundingBox(const TmplParticle &particle)
      {
        Vec3Py position = Vec3Py(particle.getPos());

        return m_bbox.intersectsWithVec3Py(position);
      }

      boost::python::list m_pyList;
      BoundingBoxPy m_bbox;
    };

    /*!
     * Class for visiting the particles to calculate their current minimum and maximum extents.
     *
     */
    class MinMaxVisitor
    {
    public:
      MinMaxVisitor(std::vector<int> dim)
        :
        m_dbl_NaN(std::numeric_limits<double>::quiet_NaN()),
        m_cur_min_pt(m_dbl_NaN,m_dbl_NaN,m_dbl_NaN),
        m_cur_max_pt(m_dbl_NaN,m_dbl_NaN,m_dbl_NaN),
        m_particle_dimensions(dim)
      {
      }

      ~MinMaxVisitor()
      {
      }

      Vec3Py getCurMinPt()
      {
        return m_cur_min_pt;
      }

      Vec3Py getCurMaxPt()
      {
        return m_cur_max_pt;
      }

      void visitParticle(const CParticle &particle)
      {
        visitSimpleParticle(particle);
      }
      void visitRotParticle(const CRotParticle &particle)
      {
        visitSimpleParticle(particle);
      }
      void visitRotParticleVi(const CRotParticleVi &particle)
      {
        visitSimpleParticle(particle);
      }
      void visitRotThermParticle(const CRotThermParticle &particle)
      {
        visitSimpleParticle(particle);
      }

    private:
      template <typename TmplParticle>
      void visitSimpleParticle(const TmplParticle &particle)
      {
        Vec3Py position = Vec3Py(particle.getPos());
        double radius = particle.getRad();

        for (int j=0, k=m_particle_dimensions[j]; j<3 && k<3; k=m_particle_dimensions[++j]) {
          if (m_cur_min_pt[k]!=m_cur_min_pt[k] || position[k]-radius<m_cur_min_pt[k])
            m_cur_min_pt[k] = position[k]-radius;
          if (m_cur_max_pt[k]!=m_cur_max_pt[k] || position[k]+radius>m_cur_max_pt[k])
            m_cur_max_pt[k] = position[k]+radius;
        }
      }

      // Variables for calculating the current particle endpoints in the geometry.
      double m_dbl_NaN;
      Vec3Py m_cur_min_pt;
      Vec3Py m_cur_max_pt;
      std::vector<int> m_particle_dimensions;
    };

    void LsmMpiPy::visitParticlesWithId(
      const boost::python::list &idList,
      boost::python::object &pyObject
    )
    {
      ParticleVisitor visitor(pyObject);
      getLatticeMaster().visitParticles(
        bpu::listToVector<int>(idList),
        visitor
      );
    }

    boost::python::list LsmMpiPy::getParticleList()
    {
      ListGatherVisitor visitor;
      getLatticeMaster().visitParticles(
        CLatticeMaster::IdVector(),
        visitor
      );
      return visitor.getList();
    }

    boost::python::list LsmMpiPy::getParticleWithIdList(
      const boost::python::list &idList
    )
    {
      ListGatherVisitor visitor;
      getLatticeMaster().visitParticles(
        bpu::listToVector<int>(idList),
        visitor
      );
      return visitor.getList();
    }

    boost::python::list LsmMpiPy::getParticlesInBBox(
      const BoundingBoxPy &bbox
    )
    {
      BBoxGatherVisitor visitor(bbox);
      getLatticeMaster().visitParticles(
        CLatticeMaster::IdVector(),
        visitor
      );
      return visitor.getList();
    }

    void LsmMpiPy::getInitMinMaxPt(Vec3Py &initMinPt, Vec3Py &initMaxPt)
    {
      getLatticeMaster().getInitMinMaxPt(initMinPt, initMaxPt);
    }

    void LsmMpiPy::getCurMinMaxPt(Vec3Py &curMinPt, Vec3Py &curMaxPt)
    {
      MinMaxVisitor visitor(getLatticeMaster().getParticleDimensions());
      getLatticeMaster().visitParticles(CLatticeMaster::IdVector(), visitor);
      curMinPt = visitor.getCurMinPt();
      curMaxPt = visitor.getCurMaxPt();
    }

     void LsmMpiPy::createBonds(
      const std::string          &groupName,
      const ParticleIdPairVector &idPairVector
    )
    {
      //getLatticeMaster().createBonds(groupName, idPairVector);
    }

    LsmMpiPy::ParticleIdPairVector
    LsmMpiPy::getBondGroupIdPairs(
        const std::string &groupName
    )
    {
      //return getLatticeMaster().getBondGroupIdPairs(groupName);
    }


    /*!
      export the interfaces to Python via boost
    */
    void exportLsm()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      setVerbosityPy(0);

      exportRunnable();
      exportCheckPointPrms();
      exportInteractionPrms();
      exportMeshBuildPrms();
      exportBondedTriMeshPrms();
      exportElasticTriMeshPrms();
      exportElasticMesh2DPrms();
      exportBondedMesh2dPrms();
      exportWallPrms();
      exportSphereBodyPrms();
      exportTriggerPrms();

      using boost::python::arg;
      boost::python::def(
        "checkMpiDimensions",
        &checkMpiDimensionsPy,
        (
          arg("numProcesses"),
          arg("mpiDimList")
        )
      );

      boost::python::def(
        "checkParticleType",
        &checkParticleTypePy,
        (
          arg("particleType")
        )
      );

      boost::python::def(
        "setVerbosity",
        &setVerbosityPy
      );

      boost::python::def(
        "setVerbosityLevel",
        &setVerbosityLevelPy
      );

      boost::python::class_<LsmMpiPy>(
        "LsmMpi",
        "Lattice Solid Model (parallelised using Message Passing Interface)"
        " container class for defining particle simulations.\n",
        boost::python::init<
          int,
          const python::list &
        >(
          (
            arg("numWorkerProcesses"),
            arg("mpiDimList")
          ),
          "Sets up MPI worker processes for use in a cartesian grid decomposition"
          " of the particle domain.\n"
          "@type numWorkerProcesses: int\n"
          "@kwarg numWorkerProcesses: The number of MPI worker processes.\n"
          "@type mpiDimList: list of 3 ints\n"
          "@kwarg mpiDimList: List of 3 elements specifying the decomposition of"
          " the simulation domain among worker processes. For example,"
          " C{[3,2,4]} indicates a rectangular grid where the M{x}-axis of the"
          " domain is divided into 3, the M{y}-axis is divided in 2 and the"
          "  M{z}-axis is divided into 4. A zero value in the list indicates"
          " that dimension discretization is chosen according to the"
          " numWorkerProcesses argument, for example,"
          " C{(numWorkerProcesses=16, mpiDimList=[2,2,0])} implies that the"
          " M{z}-axis of the will be divided into 4 parts (ie four two-by-two"
          " cells). Particles in each domain cell are handled by a MPI worker process.\n"
        )
      )
      /****fluid contents: begin****/
      .def(
        "createFluidInteraction",
        &LsmMpiPy::createFluidInteraction,
        (
          arg("cellside"),
          arg("Bw"),
          arg("Bp"),
          arg("Mu"),
          arg("alpha"),
          arg("flowrate"),
          arg("pressure"),
          arg("inflow"),
          arg("outflow"),
          arg("fluid_timestep")
        ),
        "Creates the fluid interaction. \n"
        "@type cellside: float\n"
        "@kwarg cellside: side length of fluid cells.\n"
        "@type Bw: float\n"
        "@kwarg Bw: Bulk modulus of water.\n"
        "@type Bp: float\n"
        "@kwarg Bp: Bulk modulus of particle.\n"
        "@type Mu: float\n"
        "@kwarg Mu: Viscosity of water.\n"
        "@type alpha: float\n"
        "@kwarg alpha: Adjusting factor between two time steps.\n"
        "@type flowrate: float\n"
        "@kwarg flowrate: Rate of inflow.\n"
        "@type pressure: float\n"
        "@kwarg pressure: Gradient of pressure.\n"
        "@type inflow: L{esys.lsm.util.Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg inflow: Directions of inflows.\n"
        "@type outflow: L{esys.lsm.util.Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg outflow: Directions of outflows.\n"
        "@type fluid_timestep: double\n"
        "@kwarg fluid_timestep: time step size for fluid calculations.\n"
      )
      .def(
        "createFluidInteractionVec3",
        &LsmMpiPy::createFluidInteractionVec3,
        (
          arg("cellside"),
          arg("Bw"),
          arg("Bp"),
          arg("Mu"),
          arg("alpha"),
          arg("flowrate"),
          arg("pressure"),
          arg("inflow"),
          arg("outflow"),
          arg("fluid_timestep")
        ),
        "Creates the fluid interaction. \n"
        "@type cellside: L{esys.lsm.util.Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg cellside: side length of fluid cells.\n"
        "@type Bw: float\n"
        "@kwarg Bw: Bulk modulus of water.\n"
        "@type Bp: float\n"
        "@kwarg Bp: Bulk modulus of particle.\n"
        "@type Mu: float\n"
        "@kwarg Mu: Viscosity of water.\n"
        "@type alpha: float\n"
        "@kwarg alpha: Adjusting factor between two time steps.\n"
        "@type flowrate: float\n"
        "@kwarg flowrate: Rate of inflow.\n"
        "@type pressure: float\n"
        "@kwarg pressure: Gradient of pressure.\n"
        "@type inflow: L{esys.lsm.util.Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg inflow: Directions of inflows.\n"
        "@type outflow: L{esys.lsm.util.Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg outflow: Directions of outflows.\n"
        "@type fluid_timestep: double\n"
        "@kwarg fluid_timestep: time step size for fluid calculations.\n"
      )
      .def(
        "createFieldSaver",
        &LsmMpiPy::createFluidScalarFieldSaver,
        (arg("prms"))
      )
      .def(
        "createFieldSaver",
        &LsmMpiPy::createFluidVectorFieldSaver,
        (arg("prms"))
      )
      .def(
        "createFieldSaver",
        &LsmMpiPy::createFluidInteractionScalarFieldSaver,
        (arg("prms"))
      )
      .def(
        "createFieldSaver",
        &LsmMpiPy::createFluidInteractionVectorFieldSaver,
        (arg("prms"))
      )
      /****fluid contents: end****/
      .def(
        "getNumWorkerProcesses",
        &LsmMpiPy::getNumWorkerProcesses,
        "@rtype: int\n"
        "@return: Number of spawned MPI worker processes."
      )
      .def(
        "initVerletModel",
        &LsmMpiPy::initVerletModel,
        (
          arg("particleType"),
          arg("gridSpacing"),
          arg("verletDist")
        ),
        "Initialises simulation data structures.\n"
        "@type particleType: string\n"
        "@kwarg particleType: Specifies the type discrete-element"
        " particles.\n"
        "@type gridSpacing: float\n"
        "@kwarg gridSpacing: The size of the grid-spacing used by the"
        " contact detection algorithm. For spherical particles,"
        " C{gridSpacing} must be greater than double the maximum particle radius.\n"
        "@type verletDist: float\n"
        "@kwarg verletDist: When the magnitude of a"
        " particle-displacement exceeds this amount, the neighbour lists"
        " (used in the contact detection) are updated.\n"
      )
      .def(
        "initNeighborSearch",
        &LsmMpiPy::initVerletModel,
        (
          arg("particleType"),
          arg("gridSpacing"),
          arg("verletDist")
        ),
        "Initialises simulation data structures.\n"
        "@type particleType: string\n"
        "@kwarg particleType: Specifies the type discrete-element"
        " particles.\n"
        "@type gridSpacing: float\n"
        "@kwarg gridSpacing: The size of the grid-spacing used by the"
        " contact detection algorithm. For spherical particles,"
        " C{gridSpacing} must be greater than double the maximum particle radius.\n"
        "@type verletDist: float\n"
        "@kwarg verletDist: When the magnitude of a"
        " particle-displacement exceeds this amount, the neighbour lists"
        " (used in the contact detection) are updated.\n"
      )
      .def(
        "initNeighbourSearch",
        &LsmMpiPy::initVerletModel,
        (
          arg("particleType"),
          arg("gridSpacing"),
          arg("verletDist")
        ),
        "Initialises simulation data structures.\n"
        "@type particleType: string\n"
        "@kwarg particleType: Specifies the type discrete-element"
        " particles.\n"
        "@type gridSpacing: float\n"
        "@kwarg gridSpacing: The size of the grid-spacing used by the"
        " contact detection algorithm. For spherical particles,"
        " C{gridSpacing} must be greater than double the maximum particle radius.\n"
        "@type verletDist: float\n"
        "@kwarg verletDist: When the magnitude of a"
        " particle-displacement exceeds this amount, the neighbour lists"
        " (used in the contact detection) are updated.\n"
      )
      .def(
        "setTimeStepSize",
        &LsmMpiPy::setTimeStepSize,
        (arg("dt")),
        "Sets the size of the time step used in the integration method.\n"
        "@type dt: float\n"
        "@kwarg dt: time step size used in explicit "
        " time-stepping scheme."
      )
      .def(
        "getParticleType",
        &LsmMpiPy::getParticleType,
        "@rtype: string\n"
        "@return: A string indicating the type of discrete-element"
        " particles in the model."
      )
      .def(
        "getTimeStepSize",
        &LsmMpiPy::getTimeStepSize,
        "Returns the current time-step size.\n"
        "@rtype: float\n"
        "@return: time step size for the current time-step."
      )
      .def(
        "getNumTimeSteps",
        &LsmMpiPy::getNumTimeSteps,
        "Returns the maximum number of time-steps which will"
        " be executed by the L{run} method.\n"
        "@rtype: int\n"
        "@return: Returns the maximum number of time-steps to execute."
      )
      .def(
        "setNumTimeSteps",
        &LsmMpiPy::setNumTimeSteps,
        (arg("numTimeSteps")),
        "Sets the maximum number of time-steps for the L{run}"
        " method to execute.\n"
        "@type numTimeSteps: int\n"
        "@kwarg numTimeSteps: The number of time-steps to execute in the"
        " L{run} method."
      )
      .def(
        "getTimeStep",
        &LsmMpiPy::getTimeStep,
        "Returns the current time-step number (number of time-steps"
        " executed so far).\n"
        "@rtype: int\n"
        "@return: Number of time-steps executed."
      )
      .def(
        "readGeometry",
        &LsmMpiPy::readGeometry,
        (arg("fileName")),
        "Reads particle-model setup from file.\n"
        "@type fileName: string\n"
        "@kwarg fileName: A file containing domain, particle and"
        " connection data."
      )
      .def(
        "getNumParticles",
        &LsmMpiPy::getNumParticles,
        "Returns the current number of particles in the model.\n"
        "@rtype: int\n"
        "@return: Number of particles."
      )
      .def(
        "createConnections",
        &LsmMpiPy::createConnections,
        (arg("iterable")),
        "Method for establishing bonds between particles.\n"
        "@type iterable: iterator\n"
        "@kwarg iterable: An object supporting the I{iterator protocol}"
        " (ie supports C{iter(iterable)} which forms a sequence of "
        "L{esys.lsm.geometry.TaggedIdConnection"
        "<esys.lsm.geometry.GeometryPy.TaggedIdConnection>} objects,"
        " indicating an association between particles with the specified"
        " id's."
      )
      .def("setTimingFileName",
           &LsmMpiPy::setTimingFileName,
           (arg("fileName")),
           "Method to switch on the saving of timing information and set the filename.\n"
           "@type fileName: string\n"
           "@kwarg fileName: the name of the file to which the timing data is saved\n"
      )
      .def("setSlaveTimingFileName",
           &LsmMpiPy::setSlaveTimingFileName,
           (arg("fileName")),
           "Method to switch on the saving of timing information and set the filename prefix for the slaves.\n"
           "@type fileName: string\n"
           "@kwarg fileName: the prefix of the file name to which the timing data is saved\n"
        )
       .def(
        "createParticles",
        &LsmMpiPy::createParticles,
        (arg("iterable")),
        "Creates discrete-element particles within the model.\n"
        "@type iterable: iterator\n"
        "@kwarg iterable: An object supporting the I{iterator protocol}"
        " (ie supports C{iter(iterable)} which forms a sequence of "
        "L{esys.lsm.geometry.SimpleSphere"
        "<esys.lsm.geometry.GeometryPy.SimpleSphere>} objects."
        " A model particle is created for each element "
        " object in the sequence."
      )
      .def(
        "createParticle",
        &LsmMpiPy::createParticle,
        (arg("particle")),
        "Creates a discrete-element particle within the model.\n"
        "@type particle: object\n"
        "@kwarg particle: An object from which a particle can be constructed."
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createDamping,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createLocalDamping,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createRotLocalDamping,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createABCDamping,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createGravity,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createBuoyancy,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createNRotElasticTriMeshInteractGrp,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createNRotBondedTriMeshInteractGrp,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createNRotElasticMesh2DInteractGrp,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createNRotElasticLinMeshInteractGrp,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createNRotBondedLinMeshInteractGrp,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createNRotBondInteractGrp,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createNRotShortBondInteractGrp,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createCappedNRotBondInteractGrp,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createNRotFrictionInteractGrp,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createSpringDashpotFrictionInteractGrp,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createNRotElasticInteractGrp,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createRotBondInteractGrp,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createRotThermBondInteractGrp,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createBrittleBeamInteractGrp,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createBrittleBeamSCInteractGrp,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createBrittleBeamDZCInteractGrp,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createFrictionInteractGrp,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createRotFrictionInteractGrp,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createRotThermFrictionInteractGrp,
        (arg("prms"))
      )
      .def(
          "createInteractionGroup",
          &LsmMpiPy::createVWFrictionIG,
          (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createRotElasticInteractGrp,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createRotThermElasticInteractGrp,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createHertzianElasticIG,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createHertzianViscoElasticFrictionIG,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createHertzianViscoElasticIG,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createHertzMindlinIG,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createHertzMindlinViscoIG,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createLinearDashpotIG,
        (arg("prms"))
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createNRotElasticWall,
        (
          arg("prms")
        )
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createNRotElasticSphereBody,
        (
          arg("prms")
        )
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createNRotBondedWall,
        (
          arg("prms")
        )
      )
      .def(
        "createInteractionGroup",
        &LsmMpiPy::createNRotSoftBondedWall,
        (arg("prms")),
        "Creates a group of interactions with specified properties,"
        " or a model-wall which is described by an infinite plane.\n"
        "@type prms: L{esys.lsm.InteractionPrms"
        "<esys.lsm.LsmPy.InteractionPrms>}\n"
        "@kwarg prms: An object describing the type of interaction"
        " and any parameters associated with that interaction,"
        " or parameters describing the planar wall geometry"
        " and the interactions to which particles are subjected when"
        " they encounter the wall.  For non-rotational soft bonded walls,"
        " elastic coefficients are direction dependent.\n"
      )
      .def(
          "removeInteractionGroup",
          &LsmMpiPy::removeInteractionGrp,
          (arg("name")),
          "Removes the interaction group with the specified name\n"
          "@type name: string\n"
          "@kwarg name: name of interaction group to remove\n"
      )
      .def(
        "createDamping",
        &LsmMpiPy::createDamping,
        (arg("prms")),
        "Creates viscosity within the model.\n"
        "@type prms: L{DampingPrms}\n"
        "@kwarg prms: Object describing the type of damping and"
        " the parameters associated with the damping.\n"
        "\n@status: Deprecated, use"
        " C{createInteractionGroup(DampingPrms(...))}.\n"
      )
      .def(
        "createGravity",
        &LsmMpiPy::createGravity,
        (arg("prms")),
        "Creates a gravitational body force within the model.\n"
        "@type prms: L{GravityPrms<esys.lsm.LsmPy.GravityPrms>}\n"
        "@kwarg prms: Parameters specifying gravitational acceleration.\n"
      )
      .def(
        "createBuoyancy",
        &LsmMpiPy::createBuoyancy,
        (arg("prms")),
        "Creates a buoyancy body force within the model.\n"
        "@type prms: L{BuoyancyPrms<esys.lsm.LsmPy.BuoyancyPrms>}\n"
        "@kwarg prms: Parameters specifying buoyancy force.\n"
      )
      .def(
        "createInteractionGroupTagged",
        &LsmMpiPy::createRotFrictionInteractGrpTag,
        (
          arg("prms"),
          arg("tag1"),
          arg("mask1"),
          arg("tag2"),
          arg("mask2")
        )
      )
      .def(
        "createInteractionGroupTagged",
        &LsmMpiPy::createFrictionInteractGrpTag,
        (
          arg("prms"),
          arg("tag1"),
          arg("mask1"),
          arg("tag2"),
          arg("mask2")
        )
      )
      .def(
        "createInteractionGroupTagged",
        &LsmMpiPy::createNRotFrictionInteractGrpTag,
        (
          arg("prms"),
          arg("tag1"),
          arg("mask1"),
          arg("tag2"),
          arg("mask2")
        )
      )
      .def(
        "createInteractionGroupTagged",
        &LsmMpiPy::createSpringDashpotFrictionInteractGrpTag,
        (
          arg("prms"),
          arg("tag1"),
          arg("mask1"),
          arg("tag2"),
          arg("mask2")
        )
      )
      .def(
        "createInteractionGroupTagged",
        &LsmMpiPy::createLinearDashpotInteractGrpTag,
        (
          arg("prms"),
          arg("tag1"),
          arg("mask1"),
          arg("tag2"),
          arg("mask2")
        )
      )
      .def(
        "createInteractionGroupTagged",
        &LsmMpiPy::createElasticInteractGrpTag,
        (
          arg("prms"),
          arg("tag1"),
          arg("mask1"),
          arg("tag2"),
          arg("mask2")
        )
      )
      .def(
        "createInteractionGroupTagged",
        &LsmMpiPy::createRotElasticInteractGrpTag,
        (
          arg("prms"),
          arg("tag1"),
          arg("mask1"),
          arg("tag2"),
          arg("mask2")
        ),
        "Creates a tagged interaction group between particles with the specified tags\n"
        "@type prms: L{NRotFrictionPrms}, L{NRotElasticPrms}, L{LinearDashpotPrms}, L{RotElasticPrms} etc.\n"
        "@kwarg prms: Parameters defining type of interaction to create\n"
        "@type tag1: int\n"
        "@kwarg tag1: tag of first group of particles involved in interaction\n"
        "@type mask1: int\n"
        "@kwarg mask1: mask of first group of particles involved in interaction\n"
        "@type tag2: int\n"
        "@kwarg tag2: tag of second group of particles involved in interaction\n"
        "@type mask2: int\n"
        "@kwarg mask2: mask of second group of particles involved in interaction\n"
      )
      .def(
        "createInteractionGroupTagged",
        &LsmMpiPy::createNRotElasticWallTagged,
        (
          arg("prms"),
          arg("tag"),
          arg("mask")
        ),
        "Creates an elastic interaction group between particles with the specified tag and a wall\n"
        "@kwarg prms: Parameters defining interaction to create\n"
        "@type tag: int\n"
        "@kwarg tag: tag of particles involved in interaction\n"
        "@type mask: int\n"
        "@kwarg mask: mask for tags of particles involved in interaction\n"
      )
      .def(
        "createExclusion",
        &LsmMpiPy::createExclusion,
        (
          arg("interactionName1"),
          arg("interactionName2")
        ),
        "Creates an I{interaction exclusion}. When a pair of particles"
        " come into contact, a decision is made as to the types of "
        " interactions to which the pair are subjected."
        " An exclusion precludes a pair of particles from being"
        " subjected to the interaction C{interactionName2} if"
        " they are already subjected to C{interactionName1}."
        " This is necessary, for instance, in fracture models, where"
        " particles are initially elastically bonded, but after fracture"
        " occurs particles are subjected to a frictional type of"
        " interaction.\n"
        "@type interactionName1: string\n"
        "@kwarg interactionName1: Name of an existing interaction.\n"
        "@type interactionName2: string\n"
        "@kwarg interactionName2: Name of an existing interaction.\n"
      )
      .def(
        "createSphereBody",
        &LsmMpiPy::createSphereBody,
        (
          arg("name"),
          arg("posn"),
          arg("radius")
        ),
        "Creates a sphere body -- a rigid, non-inertial sphere"
        " similar to a planar wall.\n"
        "@type name: str\n"
        "@kwarg name: Name assigned to the sphere body, used to"
        " reference this sphere body in other sphere body-related methods.\n"
        "@type posn: L{esys.lsm.util.Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg posn: The initial centre point of the sphere body.\n"
        "@type radius: double\n"
        "@kwarg radius: The radius of the sphere body.\n"
      )
      .def(
        "createWall",
        &LsmMpiPy::createWall,
        (
          arg("name"),
          arg("posn"),
          arg("normal")
        ),
        "Creates a model wall, an infinite plane specified as point-on-plane"
        " and normal-to-plane.\n"
        "@type name: str\n"
        "@kwarg name: Name assigned to the wall. This name can be used to"
        "to reference this wall in other wall-related methods.\n"
        "@type posn: L{esys.lsm.util.Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg posn: A point on the plane which describes the position of"
        " the wall.\n"
        "@type normal: L{esys.lsm.util.Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg normal: The normal to the plane describing the orientation of"
        " the wall.\n"
      )
      .def(
        "getWallPosition",
        &LsmMpiPy::getWallPosition,
        (
          arg("name")
        ),
        "Get position of named wall. Returns (0,0,0) for unknown walls.\n"
        "@type name: str\n"
        "@kwarg name: Name of the wall.\n"
      )
      .def(
        "getWallForce",
        &LsmMpiPy::getWallForce,
        (
          arg("name")
        ),
        "Get force acting on named wall. Returns (0,0,0) for unknown walls.\n"
        "@type name: str\n"
        "@kwarg name: Name of the wall.\n"
      )
      .def(
        "getSphereBodyPosition",
        &LsmMpiPy::getSphereBodyPosition,
        (
          arg("name")
        ),
        "Get position of named sphere body. Returns (0,0,0) for unknown bodies.\n"
        "@type name: str\n"
        "@kwarg name: Name of the sphere body.\n"
      )
      .def(
        "getSphereBodyForce",
        &LsmMpiPy::getSphereBodyForce,
        (
          arg("name")
        ),
        "Get force acting on named sphere body. Returns (0,0,0) for unknown bodies.\n"
        "@type name: str\n"
        "@kwarg name: Name of the sphere body.\n"
      )
      .def(
        "createRestartCheckPointer",
        &LsmMpiPy::createCheckPointer,
        (arg("prms")),
        "Causes simulation to periodically save the entire model state"
        " to file in a way neccessary to restart the simulation.\n"
        "@type prms: L{CheckPointPrms<esys.lsm.LsmPy.CheckPointPrms>}\n"
        "@kwarg prms: Object describing when and where to save model state.\n"
      )
      .def(
        "createCheckPointerWriteThroughMaster",
        &LsmMpiPy::createCheckPointerThroughMaster,
        (arg("prms")),
        "Causes simulation to periodically save the entire model state"
        " to file in a way neccessary to restart the simulation.\n"
        "@type prms: L{CheckPointPrms<esys.lsm.LsmPy.CheckPointPrms>}\n"
        "@kwarg prms: Object describing when and where to save model state.\n"
      )
      .def(
        "createCheckPointer",
        &LsmMpiPy::createSnapShots,
        (arg("prms")),
        "Causes simulation to periodically save the entire model state"
        " to file for post-processing / visualisation.\n"
        "@type prms: L{CheckPointPrms<esys.lsm.LsmPy.CheckPointPrms>}\n"
        "@kwarg prms: Object describing when and where to save model state.\n"
      )
      .def(
        "loadCheckPoint",
        &LsmMpiPy::loadCheckPoint,
        (arg("filename")),
        "Load stored model state from a checkpoint file\n"
        "@type prms: string\n"
        "@kwarg prms: name of the checkpoint file.\n"
      )
      .def(
        "readMesh",
        &LsmMpiPy::readMesh,
        (
          arg("fileName"),
          arg("meshName")
        ),
        "Creates a 3D triangulated surface mesh within the model which"
        " is loaded from file.\n"
        "@type fileName: string\n"
        "@kwarg fileName: Name of file from which surface-mesh is read.\n"
        "@type meshName: string\n"
        "@kwarg meshName: Name assigned to the created mesh.\n"
      )
      .def(
        "readMesh",
        &LsmMpiPy::readMeshWithTag,
        (
          arg("fileName"),
          arg("meshName"),
          arg("tag")
        ),
        "@type tag: int\n"
        "@kwarg tag: Only elements with tag C{tag} are created in the model"
        " (optional).\n"
      )
      .def(
        "createTriMesh",
        &LsmMpiPy::createTriMesh,
        (
          arg("meshName"),
          arg("nodeSequence"),
          arg("faceSequence")
        ),
        "Creates a 3D triangulated surface mesh within the model.\n"
        "@type meshName: string\n"
        "@kwarg meshName: Name assigned to the created mesh.\n"
        "@type nodeSequence: indexable\n"
        "@kwarg nodeSequence: Indexable sequence of node data, each element is"
        " a tuple C{(nodeId,coordinate,optionalNodeTag)}.\n"
        "@type faceSequence: indexable\n"
        "@kwarg faceSequence: Indexable sequence of triangle data, each"
        "  element is a tuple"
        " C{(faceId,(nodeId0,nodeId1,nodeId2),optionalFaceTag)}.\n"
      )
      .def(
        "readMesh2D",
        &LsmMpiPy::readMesh2D,
        (
          arg("fileName"),
          arg("meshName"),
          arg("tag")
        ),
        "Creates a 2D I{linear} surface mesh within the model which"
        " is loaded from file.\n"
        "@type fileName: string\n"
        "@kwarg fileName: Name of file from which surface-mesh is read.\n"
        "@type meshName: string\n"
        "@kwarg meshName: Name assigned to the created mesh.\n"
        "@type tag: int\n"
        "@kwarg tag: Only elements with tag C{tag} are created in the model."
      )
      .def(
        "moveSingleMeshNodeBy",
        &LsmMpiPy::moveSingleMeshNodeBy,
        (
          arg("meshName"),
          arg("nodeId"),
          arg("delta")
        ),
        "Moves an individual mesh node/vertex.\n"
        "@type meshName: string\n"
        "@kwarg meshName: Name which identifies an existing mesh.\n"
        "@type nodeId: int\n"
        "@kwarg nodeId: The identifier of the node which is to be moved.\n"
        "@type delta: L{Vec3<esys.lsm.util.FoundationPy>}\n"
        "@kwarg delta: The node is moved by this amount.\n"
      )
      .def(
        "translateMeshBy",
        &LsmMpiPy::translateMesh,
        (
          arg("meshName"),
          arg("translation")
        ),
        "Rigidly translates a whole mesh \n"
        "@type meshName: string\n"
        "@kwarg meshName: Name which identifies an existing mesh.\n"
        "@type translation: L{Vec3<esys.lsm.util.FoundationPy>}\n"
        "@kwarg translation: The mesh is translated by this vector.\n"
      )
      .def(
        "rotateMeshBy",
        &LsmMpiPy::rotateMesh,
        (
          arg("meshName"),
	  arg("origin"),
          arg("axis"),
	  arg("angle")
        ),
        "Rigidly rotates a whole mesh \n"
        "@type meshName: string\n"
        "@kwarg meshName: Name which identifies an existing mesh.\n"
        "@type origin: L{Vec3<esys.lsm.util.FoundationPy>}\n"
        "@kwarg origin: A point on the rotation axis.\n"
	"@type axis: L{Vec3<esys.lsm.util.FoundationPy>}\n"
        "@kwarg axis: The orientation of the rotation axis.\n"
	"@type angle: double\n"
        "@kwarg axis: The rotation angle in radians.\n"
      )
      .def(
        "setWallNormal",
        &LsmMpiPy::setWallNormal,
        (
          arg("wallName"),
          arg("norm")
        ),
        "Sets a wall Normal to the supplied vector.\n"
        "@type wallName: str\n"
        "@kwarg wallName: Name of the wall whose normal is to be changed.\n"
        "@type d: L{Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg d: The new normal vector of the wall."
      )
      .def(
        "getInitMinMaxPt",
        &LsmMpiPy::getInitMinMaxPt,
        (
          arg("initMinPt"),
          arg("initMaxPt")
        ),
       "Provides the initial minimum and maximum extents of all the particles read in from a geometry file.\n"
        "@type initMinPt: L{Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg initMinPt: Placeholder for the initial minimum extent of all the particles.\n"
        "@type initMaxPt: L{Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg initMaxPt: Placeholder for the initial maximum extent of all the particles.\n"
      )
      .def(
        "getCurMinMaxPt",
        &LsmMpiPy::getCurMinMaxPt,
        (
          arg("curMinPt"),
          arg("curMaxPt")
        ),
       "Provides the current minimum and maximum extents of all the particles.\n"
        "@type curMinPt: L{Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg curMinPt: Placeholder for the current minimum extent of all the particles.\n"
        "@type curMaxPt: L{Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg curMaxPt: Placeholder for the current maximum extent of all the particles.\n"
      )
      .def(
        "setSpatialDomain",
        &LsmMpiPy::setSpatialDomain,
        (
          arg("minPt"),
          arg("maxPt")
        )
      )
      .def(
        "setSpatialDomain",
        &LsmMpiPy::setBBoxSpatialDomain,
        (arg("bBox"))
      )
      .def(
        "setSpatialDomain",
        &LsmMpiPy::setBBoxSpatialDomainWithCirc,
        (
          arg("bBox"),
          arg("circDimList")
        ),
        "Defines the rectangular particle domain for this model.\n"
        "@type bBox: L{BoundingBox<esys.lsm.util.FoundationPy.BoundingBox>}\n"
        "@kwarg bBox: Defines rectangular domain.\n"
        "@type circDimList: list of 3 bool\n"
        "@kwarg circDimList: List of 3 boolean elements indicating in which"
        " dimension the circular boundary occurs. For example,"
        " C{[True,False,False]} indicates a circular boundary at the"
        " M{x=}C{bBox.getMinPt()[0]} and M{x=}bBox.getMaxPt()[0]"
        " planes of the rectangular domain specified by C{bBox}.\n"
      )
      .def(
        "findClosestParticle",
        &LsmMpiPy::findClosestParticle,
        (arg("pt")),
        "Returns the Id of the particle closest to a specified point.\n"
        "@type posn: L{Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg posn: Finds particle closest to this point.\n"
        "@rtype: int\n"
        "@return: Id of particle closest to C{posn}."
      )
      .def(
        "getParticlePosn",
        &LsmMpiPy::getParticlePosn,
        (arg("id")),
        "Returns the location of a particle with a specified Id.\n"
        "@type id: int\n"
        "@kwarg id: Particle Id.\n"
        "@rtype: L{Vec3<esys.lsm.util.FoundationPy>}\n"
        "@return: Position of particle."
      )
      .def(
        "changeRadiusBy",
        &LsmMpiPy::changeRadiusBy,
        (
          arg("tag"),
          arg("deltaR")
        ),
        "Changes the radius of particles with a given tag.\n"
        "@type tag: int\n"
        "@kwarg tag: tag of the particle whose radius will change.\n"
        "@type deltaR: double\n"
        "@kwarg deltaR: amount by which to change the radius of particles."
      )
      .def(
        "moveParticleTo",
        &LsmMpiPy::moveSingleParticleTo,
        (
          arg("id"),
          arg("posn")
        ),
        "Sets the absolute position of a particle.\n"
        "@type id: int\n"
        "@kwarg id: Id of the particle to be moved.\n"
        "@type posn: L{Vec3<esys.lsm.util.FoundationPy>}\n"
        "@kwarg posn: New position of particle."
      )
      .def(
        "moveTaggedParticlesBy",
        &LsmMpiPy::moveTaggedParticlesBy,
        (
          arg("tag"),
          arg("displacement")
        ),
        "Translates the position of a particle by a specified displacment.\n"
        "@type tag: int\n"
        "@kwarg tag: Tag of all particles which are to be moved.\n"
        "@type displacement: L{Vec3<esys.lsm.util.FoundationPy>}\n"
        "@kwarg displacement: Amount by which particle positions are translated."
      )
      .def(
        "moveSphereBodyBy",
        &LsmMpiPy::moveSphereBodyBy,
        (
          arg("sphereName"),
          arg("d")
        ),
        "Moves a sphere body by a specified displacement.\n"
        "@type sphereName: str\n"
        "@kwarg sphereName: Name of the sphere body which is to be moved.\n"
        "@type d: L{Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg d: The displacement vector by which the sphere body\n"
        " is displaced. \n"
      )
      .def(
        "moveWallBy",
        &LsmMpiPy::moveWallBy,
        (
          arg("wallName"),
          arg("d")
        ),
        "Moves a model wall by a specified displacement.\n"
        "@type wallName: str\n"
        "@kwarg wallName: Name of the wall which is to be moved.\n"
        "@type d: L{Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg d: The displacement vector by which the wall is displaced.\n"
      )
      .def(
        "applyForceToWall",
        &LsmMpiPy::applyForceToWall,
        (
          arg("interactionName"),
          arg("force")
        ),
        "A wall is displaced in the normal direction so that the total"
        " force applied to the wall attains a specified value.\n"
        "@type interactionName: str\n"
        "@kwarg interactionName: Name of the interaction group governing"
        " interactions between the wall and particles.\n"
        "@type force: L{Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg force: Wall is displaced in the direction of the normal,"
        " so that net wall force is C{force.n}."
      )
      .def(
        "setParticleVelocity",
              &LsmMpiPy::setParticleVel,
        (
          arg("id"),
          arg("Velocity")
        ),
        "Set the velocity of a particle \n"
        "@type id: int\n"
        "@kwarg id: the ID of the particle \n"
        "@type Velocity: L{Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg Velocity: The velocity of the particle.\n"
      )
      .def(
        "setParticleAngularVelocity",
        &LsmMpiPy::setParticleAngVel,
        (
          arg("id"),
          arg("angularVelocity")
        ),
        "Set the angular velocity of a rotational particle \n"
        "@type id: int\n"
        "@kwarg id: the ID of the particle \n"
        "@type angularVelocity: L{Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg angularVelocity: The angular velocity of the particle.\n"
      )
      .def(
        "setParticleDensity",
        &LsmMpiPy::setParticleDensity,
        (
          arg("tag"),
          arg("mask"),
          arg("Density")
        ),
        "Set the density of a group of tagged particles \n"
        "@type tag: int\n"
        "@kwarg tag: the tag of the particles \n"
        "@type mask: int\n"
        "@kwarg mask: the tag mask \n"
        "@type Density: float\n"
        "@kwarg Density: The density of the particle.\n"
      )
     .def(
        "resetParticleOrientation",
        &LsmMpiPy::resetParticleOrientation,
        (
          arg("tag"),
          arg("mask")
        ),
        "Reset the orientation of a group of tagged particles \n"
        "@type tag: int\n"
        "@kwarg tag: the tag of the particles \n"
        "@type mask: int\n"
        "@kwarg mask: the tag mask \n"
      )
     .def(
        "setTaggedParticleVelocity",
        &LsmMpiPy::setTaggedParticleVel,
        (
          arg("tag"),
          arg("Velocity")
        ),
        "Set the velocity of a group of tagged particles \n"
        "@type tag: int\n"
        "@kwarg tag: the tag of the particles \n"
        "@type Velocity: L{Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg Velocity: The velocity of the particle.\n"
      )
      .def(
        "setVelocityOfWall",
        &LsmMpiPy::setVelocityOfWall,
        (arg("name"),arg("velocity")),
        "Set the velocity of a wall with viscous drag. This does not"
        " influence the position of the wall, only the viscous drag"
        " applied to particles interacting with the wall. Therefore"
        " it is meaningless for walls without viscous drag.\n"
        "@type name:  str\n"
        "@kwarg name: the name of the wall\n"
        "@type velocity: L{Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg velocity: The velocity of the wall.\n"
      )
      .def(
        "tagParticleNearestTo",
        &LsmMpiPy::tagParticleNearestTo,
        (arg("tag"),arg("mask"),arg("Position")),
    "Finds the particle nearest a specified position and sets its tag to the value provided\n"
    "@type tag: int\n"
    "@kwarg tag: the tag to assign to the particle\n"
    "@type mask: int\n"
    "@kwarg mask: the tag mask to apply when tagging the particle\n"
    "@type Position: C{vec3}\n"
    "@kwarg Position: the position nearest the particle to tag\n"
      )
      .def(
        "setParticleNonDynamic",
        &LsmMpiPy::setParticleNonDynamic,
        (arg("tag")),
        "Make the particle non-dynamic, i.e. the particle still interacts \n"
        "with other particles the usual way but doesn't move in response to \n"
        "forces applied to it. Useful if the particle is moved in order to \n"
        "create a deformation source.\n"
        "@type tag: int\n"
        "@kwarg tag: the tag of the particle.\n"
       )
       .def(
        "setParticleNonRotational",
        &LsmMpiPy::setParticleNonRot,
        (arg("tag")),
        "Set the particle with tag C{tag} to be non-rotational,"
        " i.e. it still participates in \n"
        "rotational (RotFriction, RotBonded, etc) interactions but doesn't"
        "  rotate in response to applied torque. Only applicable if the"
        "  particle type is rotational.\n"
        "@type tag: int\n"
        "@kwarg tag: the tag of the particle."
      )
       .def(
        "setParticleNonTranslational",
        &LsmMpiPy::setParticleNonTrans,
        (arg("tag")),
        "Set the particle with tag C{tag} to be non-translational,"
        " i.e. it still has rotational degrees of freedom\n"
        "but doesn't translate in response to forces. Only applicable if the"
        "  particle type is rotational.\n"
        "@type tag: int\n"
        "@kwarg tag: the tag of the particle."
      )
      .def(
        "addPreTimeStepRunnable",
        &LsmMpiPy::addPreTimeStepRunnable,
        (arg("runnable")),
        "Adds a L{Runnable} object to the list of runnables. The L{Runnable.run}\n"
        "method is called before the execution of a time-step. This method can be\n"
        "used to introduce a loading mechanism (apply force to walls,"
        " move walls, etc).\n"
        "@type runnable: L{Runnable<esys.lsm.LsmPy.Runnable>}\n"
        "@kwarg runnable: An object which has a C{run} method which takes"
        " no arguments.\n"
      )
      .def(
        "addPostTimeStepRunnable",
        &LsmMpiPy::addPostTimeStepRunnable,
        (arg("runnable")),
        "Like the L{addPreTimeStepRunnable} method, except the C{run} method"
        " of C{runnable} gets called at the end of a time-step.\n"
        "@type runnable: L{Runnable<esys.lsm.LsmPy.Runnable>}\n"
        "@kwarg runnable: An object which has a C{run} method which takes"
        " no arguments.\n"
      )
      .def(
        "createFieldSaver",
        &LsmMpiPy::createParticleScalarFieldSaver,
        (arg("prms"))
      )
      .def(
        "createFieldSaver",
        &LsmMpiPy::createParticleVectorFieldSaver,
        (arg("prms"))
      )
      .def(
        "createFieldSaver",
        &LsmMpiPy::createInteractionScalarFieldSaver,
        (arg("prms"))
      )
      .def(
        "createFieldSaver",
        &LsmMpiPy::createCheckedInteractionScalarFieldSaver,
        (arg("prms"))
      )
      .def(
        "createFieldSaver",
        &LsmMpiPy::createInteractionVectorFieldSaver,
        (arg("prms"))
      )
     .def(
        "createFieldSaver",
        &LsmMpiPy::createCheckedInteractionVectorFieldSaver,
        (arg("prms"))
      )
      .def(
        "createFieldSaver",
        &LsmMpiPy::createTaggedParticleScalarFieldSaver,
        (arg("prms"))
      )
      .def(
        "createFieldSaver",
        &LsmMpiPy::createTaggedParticleVectorFieldSaver,
        (arg("prms"))
      )
      .def(
        "createFieldSaver",
        &LsmMpiPy::createTaggedInteractionScalarFieldSaver,
        (arg("prms"))
      )
      .def(
        "createFieldSaver",
        &LsmMpiPy::createParticleVectorFieldSaverWithTrigger,
        (arg("prms"))
      )
      .def(
        "createFieldSaver",
        &LsmMpiPy::createTaggedParticleVectorFieldSaverWithTrigger,
        (arg("prms"))
      )
      .def(
        "createFieldSaver",
        &LsmMpiPy::addVectorWallField,
        (arg("prms")),
        "Causes specified data to be saved periodically to file.\n"
        "@type prms: L{FieldSaverPrms}\n"
        "@kwarg prms: An object specifing where, when and what data are"
        " to be saved.\n"
      )
#if 0

        .def(
          "addTaggedScalarParticleDistributionSaver",
          &LsmMpiPy::addTaggedScalarParticleDistributionSaver
        )
#endif
      .def(
        "addVectorTriangleSaveField",
        &LsmMpiPy::addVectorTriangleSaveField,
        (arg("prms")),
        "Add a vector triangle field saver to the simulation\n"
        "@type prms: L{TriangleVectorFieldSaverPrms}\n"
        "@kwarg prms: An object specifing where, when and what data are"
        " to be saved.\n"
      )
      .def(
        "addScalarTriangleSaveField",
        &LsmMpiPy::addScalarTriangleSaveField,
        (arg("prms")),
        "Add a scalar triangle field saver to the simulation\n"
        "@type prms: L{TriangleScalarFieldSaverPrms}\n"
        "@kwarg prms: An object specifing where, when and what data are"
        " to be saved.\n"
      )
      .def(
        "force2dComputations",
        &LsmMpiPy::force2dComputations,
        "Ensures particles only move in the M{x-y} plane and that"
        " rotations only occur about the M{z}-axis.\n"
      )
      .def(
        "runTimeStep",
        &LsmMpiPy::runTimeStep,
        "Runs a single step of the time-integration method.\n"
      )
      .def(
        "run",
        &LsmMpiPy::run,
        "Runs multiple steps of the time-integration method.\n"
      )
      .def(
        "exit",
        &LsmMpiPy::exit,
        "Exits the simulation after running a series of single"
        "steps of the time-integration method.\n"
      )
      .def(
        "visitNodeRefs2d",
        &LsmMpiPy::visitNodeRefs2d,
        (
          arg("meshName"),
          arg("nodeRefVisitor")
        ),
        "Method for visiting mesh node Id's.\n"
        "@type meshName: string\n"
        "@kwarg meshName: name of the mesh whose nodes will be visited.\n"
        "@type nodeRefVisitor: L{object}\n"
        "@kwarg nodeRefVisitor: object which has a method named "
        "C{visitNodeRef} which takes a single integer argument.\n"
      )
      .def(
        "visitRefStressPairs2d",
        &LsmMpiPy::visitRefStressPairs2d,
        (
          arg("meshName"),
          arg("refStressVisitor")
        ),
        "Method for visiting mesh C{(element-Id, stress)} pairs.\n"
        "@type meshName: string\n"
        "@kwarg meshName: name of the mesh whose elements will be visited.\n"
        "@type refStressVisitor: L{object}\n"
        "@kwarg refStressVisitor: object which has a method named "
        "C{visitRefStressPair} which two arguments, an integer element-id"
        " argument and a L{Vec3<esys.lsm.util.FoundationPy>} stress"
        " argument.\n"
      )
      .def(
        "visitParticlesWithId",
        &LsmMpiPy::visitParticlesWithId,
        (
          arg("idList"),
          arg("particleVisitor")
        ),
        "Method for visiting particle data.\n"
        "@type idList: list\n"
        "@kwarg idList: List of particle-id values specifying which particles"
        " are to be visited.\n"
        "@type particleVisitor: L{object}\n"
        "@kwarg particleVisitor: object which has a method named "
        "C{visitParticle} which accepts a single particle argument."
      )
      .def(
        "getParticleList",
        &LsmMpiPy::getParticleList
      )
      .def(
        "getParticleList",
        &LsmMpiPy::getParticleWithIdList,
        (
          arg("idList")
        ),
        "Returns a I{copy} of particles with specified id's. If the C{idList}"
        " argument is not specified (or is empty) then all particles are"
        " returned in the list.\n"
        "@type idList: list\n"
        "@kwarg idList: Optional list of particle-id values specifying which"
        " particles are to be returned.\n"
        "@rtype: list\n"
        "@return: Python C{list} containing particles with specified id."
        " If C{idList} is empty, then all particles are returned in the list.\n"
      )
      .def(
        "getParticlesInBBox",
        &LsmMpiPy::getParticlesInBBox,
        (
          arg("bbox")
        ),
        "Returns a I{copy} of particles currently residing within the specified"
        " C{BoundingBoxPy}.\n"
        "@type bbox: C{BoundingBoxPy}\n"
        "@kwarg bbox: Bounding Box of subregion containing the"
        " particles to be returned.\n"
        "@rtype: list\n"
        "@return: Python C{list} containing particles within"
        " specified bounding box.\n"
      )

      .def(
        "getLsmVersion",
        &LsmMpiPy::getLsmVersion,
        "Returns " PACKAGE_NAME " version string.\n"
        "@rtype: string\n"
        "@return: " PACKAGE_NAME " version."
      )
      .def(
        "setVerbosity",
        &LsmMpiPy::SetVerbosityPy,
        (
          arg("verbose")
        ),
        "Set verbosity for Master and Workers\n"
        "@type verbose: bool\n"
        "@kwarg verbose: set the amount of debug info written to file false->nothing, true->all \n"
        " For a more fine grained control use setVerbosityLevelPy which takes an int argument (0-7) \n"
      )
      .def(
        "setConsoleFilename",
        &LsmMpiPy::SetConsoleFilenamePy,
        (
          arg("filename")
        ),
        "Set console output file name for Master and Workers\n"
        "@type filename: string\n"
        "@kwarg filename: set basic file name of the console output \n"
        " The full filename for each process will be filename.rank \n"
      )
      .def(
        "setConsoleBuffered",
        &LsmMpiPy::SetConsoleBufferedPy,
        (
          arg("bsize")
        ),
        "Set console buffer mode and size for Master and Workers\n"
        "@type bsize: int\n"
        "@kwarg bsize: set the size of the console buffer \n"
        " A larger buffer may improve performance but risks losing output in case of program termination\n"
      )
      .def(
        "setVerbosityLevel",
        &LsmMpiPy::SetVerbosityLevelPy,
        (
          arg("verbose")
        ),
        "Set verbosity for Master and Workers\n"
        "@type verbose: int\n"
        "@kwarg verbose: set the amount of debug info written to file 0->nothing, 7->all \n"
      )
      .def(
        "setInteractionParameter",
        &LsmMpiPy::setInteractionParameter,
        (
          arg("igname"),
          arg("pname"),
          arg("val")
          
        ),
        "Modify a parameter of an existing interaction group during a simulation\n"
        "@type igname: string\n"
        "@kwarg igname: the name of the interaction group \n" 
        "@type pname: string\n"
        "@kwarg pname: the name of the parameter \n" 
        "@type val: float\n"
        "@kwarg val: the value by which the parameter is modified \n" 
      )
      ;
    }
  }
}
