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

#ifndef __LATTICEMASTER_H
#define __LATTICEMASTER_H

//--- Project includes ---

#include "Parallel/mpibuf.h"
#include "Parallel/mpivbuf.h"
#include "Parallel/LatticeParam.h"
#include "Parallel/RankAndComm.h"
#include "Parallel/sublattice_cmd.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "Foundation/console.h"
#include "Foundation/Runnable.h"

#include "Fields/FieldMaster.h"
#include "Fields/MaxTrigger.h"

#include "Geometry/GeometryInfo.h"

using esys::lsm::GeometryInfo;

#include "Model/Damping.h"
#include "Model/LocalDamping.h"
#include "Model/ABCDampingIGP.h"
#include "Model/Particle.h"
#include "Model/RotParticle.h"
#include "Model/RotParticleVi.h"
#include "Model/RotThermParticle.h"
#include "Model/ElasticInteraction.h"
#include "Model/FrictionInteraction.h"
#include "Model/SpringDashpotFrictionInteraction.h"
#include "Model/FractalFriction.h"
#include "Model/AdhesiveFriction.h"
#include "Model/HertzianElasticInteraction.h"
#include "Model/HertzianViscoElasticFrictionInteraction.h"
#include "Model/HertzianViscoElasticInteraction.h"
#include "Model/HertzMindlinInteraction.h"
#include "Model/HertzMindlinViscoInteraction.h"
#include "Model/LinearDashpotInteraction.h"
#include "Model/MeshData.h"
#include "Model/ETriMeshIP.h"
#include "Model/BTriMeshIP.h"
#include "Model/BMesh2DIP.h"
#include "Model/BondedInteraction.h"
#include "Model/CappedBondedInteraction.h"
#include "Model/RotFricInteraction.h"
#include "Model/RotBondedInteraction.h"
#include "Model/BrittleBeamSC.h"
#include "Model/BrittleBeamDZC.h"
#include "Model/RotElasticInteraction.h"
#include "Model/RotThermFricInteraction.h"
#include "Model/RotThermBondedInteraction.h"
#include "Model/RotThermElasticInteraction.h"
#include "Model/BodyForceGroup.h"
#include "Model/EWallInteractionGroup.h"
#include "Model/BWallInteractionGroup.h"
#include "Model/ViscWallIG.h"
#include "Model/SoftBWallInteractionGroup.h"
#include "Model/ESphereBodyInteractionGroup.h"
#include "Fields/FluidFieldMaster.h" //fluid contents
#include "Fields/FluidInteractionFieldMaster.h" //fluid contents
#include "Parallel/GMRESSolverMaster.h" //fluid contents

#include <boost/filesystem/path.hpp>

//--- MPI includes ---
#include <mpi.h>

//--- TML includes ---
#include "tml/comm/comm_world.h"

// -- STL includes --
#include <vector>
#include <list>
#include <map>
#include <utility>
#include <string>
#include <limits>

// forward decls.
// includes are in the .cpp
class CheckPointController;

namespace esys
{
  namespace lsm
  {
    class GeometryInfo;
    class BodyForceIGP;
  }
}

/*!
  Class for the initialisation and control of a group of sublattices

  \author Steffen Abe, David Place, Shane Latham
  $Revision$
  $Date$
*/
class MpiWTimers;

namespace esys
{
  namespace lsm
  {
    typedef std::vector<int> IntVector;
  }
}

bool sortOnZ(const pair<Vec3,double>, const pair<Vec3,double>); //fluid contents
bool sortOnY(const pair<Vec3,double>, const pair<Vec3,double>); //fluid contents
bool sortOnX(const pair<Vec3,double>, const pair<Vec3,double>); //fluid contents

class CLatticeMaster
{
 public:
    typedef std::vector<esys::lsm::Runnable *> RunnableVector;
    typedef std::pair<int, int>                ParticleIdPair;
    typedef std::vector<ParticleIdPair>        ParticleIdPairVector;
    typedef std::vector<MeshNodeData> MeshNodeDataVector;
    typedef std::vector<MeshTriData> MeshTriDataVector;
    typedef std::pair<MeshNodeDataVector,MeshTriDataVector> TriMeshDataPair;

 private:
    std::string                           m_timingFileName;
    MpiWTimers                            *m_pTimers;
    CheckPointController                  *m_pCheckPointController; // for restart checkpoints
    CheckPointController                  *m_pSnapShotController; // for viz/analysis dumps
    esys::lsm::CLatticeParam::ProcessDims m_processDims;
    GMRESSolverMaster                     *m_master_solver; //fluid contents

 protected:
    typedef std::vector<int> ConnIdVector;
    map<int,ConnIdVector> m_temp_conn;
    vector<AFieldMaster*> m_save_fields;

    // -- variables for global model geometry
    GeometryInfo m_geo_info;
    bool m_bbx_has_been_set;
    bool m_geometry_is_initialized;
    // ----

    int m_global_rank;
    int m_global_size;
    int m_max_ts;
    int m_center_id;
    double m_total_time;
    int m_t ;
    int m_t_f ; //fluid contents
    int m_t_old ; //fluid contents
    double m_dt;
    double m_dt_f; //fluid contents
    bool m_isInitialized ;
    bool m_first_time;
    std::string m_particle_type;
    bool m_fluidinitiated; //fluid contents


    RunnableVector m_preRunnableVector;
    RunnableVector m_postRunnableVector;

    TML_Comm m_tml_global_comm;
    MPI_Comm m_global_comm, m_local_comm;
    MPI_Group m_mpi_local_group;  // needs to be member in order to free at desctruction

    // Variables for calculating the initial particle endpoints in the geometry.
    double m_dbl_NaN;
    Vec3 m_init_min_pt;
    Vec3 m_init_max_pt;
    esys::lsm::IntVector m_particle_dimensions;

    void runRunnables(RunnableVector::iterator begin, RunnableVector::iterator end);
    void runPreRunnables();
    void runPostRunnables();

    void saveTimingData();
    TriMeshDataPair readTriMesh(const std::string &fileName,int);
    TriMeshDataPair readTriMesh(const std::string &fileName);
    void readAndDistributeMesh2D(const std::string&,int);

  MpiRankAndComm getGlobalRankAndComm() const
  {
    return MpiRankAndComm(m_global_rank, m_global_comm);
  }

public:
    CLatticeMaster();
    ~CLatticeMaster();

    std::string getLsmVersion() const
    {
      return std::string(PACKAGE_VERSION);
    }

    /****fluid contents: begin****/
    void addFluidInteraction(double,double, double, double, double, double, double, Vec3, Vec3, double);
    void addFluidInteractionVec3(Vec3,double, double, double, double, double, double, Vec3, Vec3, double);
    void addScalarFluidSaveField(const string&,const string&,const string&,int,int,int);
    void addVectorFluidSaveField(const string&,const string&,const string&,int,int,int);
    void addScalarFluidInteractionSaveField(const string&,const string&,const string&,int,int,int);
    void addVectorFluidInteractionSaveField(const string&,const string&,const string&,int,int,int);
    void calculatePressure();
    vector<pair<Vec3,double> > collectSingleCoeffi(const string&);
    vector<pair<Vec3,double> > sortCoefficient(vector<pair<Vec3,double> > coefficient); 
    void distributePressure(vector<pair<Vec3,double> >);
    /****fluid contents: end****/

    int getNumWorkerProcesses() const;

    int getTimeStep() const {return m_t;}
    double getTimeStepSize() const {return m_dt;}
    void setTimeStepSize(double dt);

    void init();

    void setupWorkers(int numWorkers);
    void run();
    void runInit();
    void runOneStep();
    void runEnd();
    void oneStep();
    void searchNeighbors(bool);
    bool checkNeighbors();
    void updateInteractions();
    void addBondedIG(const CBondedIGP&);
    void addCappedBondedIG(int,const std::string&,double,double,double);
    void addShortBondedIG(int,const std::string&,double,double);

    void addPairIG(const CElasticIGP &prms);
    void addPairIG(const CFrictionIGP &prms);
    void addPairIG(const CSpringDashpotFrictionIGP &prms);
    void addPairIG(const FractalFrictionIGP &prms);
    void addPairIG(const CAdhesiveFrictionIGP &prms);
    void addPairIG(const CRotElasticIGP &prms);
    void addPairIG(const CRotFrictionIGP &prms);
    void addPairIG(const CHertzianElasticIGP &prms);
    void addPairIG(const CHertzianViscoElasticFrictionIGP &prms);
    void addPairIG(const CHertzianViscoElasticIGP &prms);
    void addPairIG(const CHertzMindlinIGP &prms);
    void addPairIG(const CHertzMindlinViscoIGP &prms);
    void addPairIG(const CLinearDashpotIGP &prms);
    void addPairIG(const CRotThermElasticIGP &prms);
    void addPairIG(const CRotThermFrictionIGP &prms);
    void addTaggedPairIG(const CRotFrictionIGP &prms,int,int,int,int);
    void addTaggedPairIG(const CFrictionIGP &prms,int,int,int,int);
    void addTaggedPairIG(const CSpringDashpotFrictionIGP &prms,int,int,int,int);
    void addTaggedPairIG(const CHertzianElasticIGP &prms,int,int,int,int);
    void addTaggedPairIG(const CHertzianViscoElasticFrictionIGP &prms,int,int,int,int);
    void addTaggedPairIG(const CHertzianViscoElasticIGP &prms,int,int,int,int);
    void addTaggedPairIG(const CHertzMindlinIGP &prms,int,int,int,int);
    void addTaggedPairIG(const CHertzMindlinViscoIGP &prms,int,int,int,int);
    void addTaggedPairIG(const CLinearDashpotIGP &prms,int,int,int,int);
    void addTaggedPairIG(const CRotElasticIGP &prms,int,int,int,int);
    void addTaggedPairIG(const CElasticIGP &prms,int,int,int,int);

    void removeIG(const std::string&);

    void readAndDistributeTriMesh(const std::string&,const std::string&,int);
    void readAndDistributeTriMesh(const std::string&,const std::string&);
    void createTriMesh(
      const std::string &meshName,
      const MeshNodeDataVector &mndVector,
      const MeshTriDataVector &mtdVector
    );
    void addMesh2D(const std::string&,const std::string&,int);
    void addMesh2DIG(const ETriMeshIP &prms);
    void addTriMesh(const std::string &meshName, const std::string &fileName);
    void addTriMeshIG(const ETriMeshIP &prms);

    void addBondedTriMeshIG(const BTriMeshIP &triMeshPrms, const MeshTagBuildPrms &buildPrms);
    void addBondedTriMeshIG(const BTriMeshIP &triMeshPrms, const MeshGapBuildPrms &buildPrms);

    void addBondedMesh2DIG(const BMesh2DIP&, const MeshTagBuildPrms&);
    void addBondedMesh2DIG(const BMesh2DIP&, const MeshGapBuildPrms&);
    void addDamping(const CDampingIGP &dampingIGP);
    void addDamping(const CLocalDampingIGP &dampingIGP);
    void addDamping(const ABCDampingIGP &dampingIGP);

    void addSingleIG(const esys::lsm::GravityIGP &gravityIGP);
    void addSingleIG(const esys::lsm::BuoyancyIGP &buoyancyIGP);
    void addExIG(const std::string&,const std::string&);
    void setNumSteps(int s);
    int getNumSteps() const {return m_max_ts;};
    int getSteps() const {return m_t;};

    void addRotBondedIG(int,const std::string&,double,double,double,double,double,double,double,double,bool,bool,double);
    void addBrittleBeamSCIG(const BrittleBeamSCIGP&);
    void addBrittleBeamDZCIG(const BrittleBeamDZCIGP&);
    void addRotThermBondedIG(const CRotThermBondedIGP &prms);

//    ParticleIdPairVector getBondGroupIdPairs(const std::string &groupName);

    // --- wall related fucntions ---
    void addWall(const std::string&,const Vec3&,const Vec3&);
    void addWallIG(const CEWallIGP&);
    void addWallIG(const CBWallIGP&);
    void addWallIG(const CVWallIGP&);
    void addWallIG(const CSoftBWallIGP&);
    void addTaggedWallIG(const CEWallIGP&,int,int);
    Vec3 getWallPosn(const std::string&);
    Vec3 getWallForce(const std::string&);

    // --- sphere body related functions ---
    void addSphereBody(const std::string&,const Vec3&,const double&);
    void addSphereBodyIG(const CESphereBodyIGP&);
    Vec3 getSphereBodyPosn(const std::string&);
    Vec3 getSphereBodyForce(const std::string&);

    //     void initSoftBondedWall(const string&,const Vec3&,const Vec3&,double,double,double,int);


    void changeRadiusBy(int particleTag, const double &deltaR);
    void moveParticleTo(int particleTag, const Vec3 &posn);
    void moveTaggedParticlesBy(int particleTag, const Vec3 &displacement);
    void moveSingleParticleTo(int particleId, const Vec3 &posn);
    Vec3 getParticlePosn(int particleId);
    void setParticleNonDynamic(int);
    void setParticleNonRot(int);
    void setParticleNonTrans(int);
    void setParticleVel(int,const Vec3&);
    void setParticleAngVel(int,const Vec3&);
    void setParticleDensity(int tag,int mask,double rho);
    void resetParticleRotation(int,int);
    void setTaggedParticleVel(int tag,const Vec3&);
    void moveWallBy(const std::string&,const Vec3&);
    void moveSphereBodyBy(const std::string&,const Vec3&);
    void setWallNormal(const std::string&,const Vec3&);
    void setVelocityOfWall(const std::string&,const Vec3&);
    void tagParticleNearestTo(int,int,const Vec3&);
    /**
     * Returns the id of the particle which is closest to the
     * specified point.
     */
    int findParticleNearestTo(const Vec3& pos);
    void applyForceToWall(const std::string&,const Vec3&);
    void applyForceToSphereBody(const std::string&,const Vec3&);
    // --- Mesh movement functions ---
    void moveSingleNodeBy(const std::string&,int,const Vec3&);
    void moveTaggedNodesBy(const std::string&,int,const Vec3&);
    void translateMeshBy(const std::string&,const Vec3&);
    void rotateMeshBy(const std::string&,const Vec3&,const Vec3&,const double);
    void applyPressureToMeshInteraction(const std::string&,double);
	

    void saveTimingDataToFile(const std::string &fileNamePrefix);

    /**
     * Enforces particles to be restricted to motion in the x-y plane.
     */
    void do2dCalculations(bool do2d);

    /**
     * Sets explicit process partitioning info for MPI_Dims_create.
     *
     * @param dims Dimension process partition vector.
     *             dims[0] x-partitioning,
     *             dims[1] y-partitioning and
     *             dims[2] z-partitioning. A zero value for a dimension
     *             is converted to a non-zero value using MPI_Dims_create
     */
    void setProcessDims(const esys::lsm::CLatticeParam::ProcessDims &dims);

    const esys::lsm::CLatticeParam::ProcessDims &getProcessDims() const;

    /*!
     * Specifies a file in which timing results are written.
     *
     * \param fileName Name of file, overwritten if it exists.
     */
    void setTimingFileName(const std::string &fileName);

    /*!
     * Returns the name of the file to which timing results are written.
     *
     * \return Name of file, empty std::string if no timing results are being
     * recorded.
     */
    const std::string &getTimingFileName() const;

    const std::string &getParticleType() const
    {
      return m_particle_type;
    }

    int getNumParticles();

    //! field saving functions
    void addScalarParticleSaveField(const std::string&,const std::string&,const std::string&,int,int,int);
    void addTaggedScalarParticleSaveField(const std::string&,const std::string&,const std::string&,int,int,int,int,int);
    void addVectorParticleSaveField(const std::string&,const std::string&,const std::string&,int,int,int);
    void addTaggedVectorParticleSaveField(const std::string&,const std::string&,const std::string&,int,int,int,int,int);
    void addScalarInteractionSaveField(const std::string&,const std::string&,const std::string&,const std::string&,const std::string&,int,int,int,bool checked=false);
    void addScalarHistoryInteractionSaveField(const std::string&,const std::string&,const std::string&,const std::string&,const std::string&,int,int,int);
    void addVectorInteractionSaveField(const std::string&,const std::string&,const std::string&,const std::string&,const std::string&,int,int,int,bool checked=false);
    void addTaggedScalarInteractionSaveField(const std::string&,const std::string&,const std::string&,const std::string&,const std::string&,int,int,int,int,int,bool);
    void addTaggedScalarParticleDistributionSaver(const std::string&,const std::string&,const std::string&,int,int,int,int,int,int,double,double,int);
    void addVectorTriangleSaveField(const string&,const string&,const string&,const string&,int,int,int);
    void addScalarTriangleSaveField(const string&,const string&,const string&,const string&,int,int,int);
    void addVectorWallField(const string&,const string&,vector<string>,const string&,int,int,int);
    // fields with trigger
    void addVectorParticleSaveFieldWT(const std::string&,const std::string&,const std::string&,int,int,int,const MaxTrigParams&);
    void addTaggedVectorParticleSaveFieldWT(const std::string&,const std::string&,const std::string&,int,int,int,int,int,const MaxTrigParams&);

    /**
     * Initialises parameters for performing model checkpointing.
     *
     * @param fileNamePrefix Path prefix for checkpoint files. Multiple checkpoint
     *                       files may be generated for a single timestep snapshot.
     * @param beginTime      Time to begin checkpointing. Time of first checkpoint.
     * @param endTime        End time for checkpointing. Time of last checkpoint.
     * @param timeInterval   Time interval between checkpoint file generation.
     * @param precision      the output precision (digits)
     */
    void performCheckPoints(
        const std::string &fileNamePrefix,
        int beginTime,
        int endTime,
        int timeInterval,
	int precision
    );

    /**
     * Initialises parameters for performing model checkpointing with writing done
     * through master process.
     *
     * @param fileNamePrefix Path prefix for checkpoint files. Multiple checkpoint
     *                       files may be generated for a single timestep snapshot.
     * @param beginTime      Time to begin checkpointing. Time of first checkpoint.
     * @param endTime        End time for checkpointing. Time of last checkpoint.
     * @param precision      the output precision (digits)
     * @param timeInterval   Time interval between checkpoint file generation.
     */

    void performCheckPointsThroughMaster(
        const std::string &fileNamePrefix,
        int beginTime,
        int endTime,
        int timeInterval,
	int precision
    );

    void initSnapShotController(const std::string&,int,int,int);

    //! initialization functions
    void makeLattice(
      const char *particleType,
      double gridSize,
      double verletDist
    );

    void makeLattice(
      const char *particleType,
      double gridSize,
      double verletDist,
      double dt
    );

    /**
     * Adds reference to an object whose 'run' method is executed at
     * the beginning of the runOneStep method.
     */
    void addPreTimeStepRunnable(esys::lsm::Runnable &runnable);

    /**
     * Returns vector of Runnable objects (which are run pre-time-step).
     */
    const RunnableVector &getPreTimeStepRunnableVector() const
    {
      return m_preRunnableVector;
    }

    /**
     * Returns vector of Runnable objects (which are run pre-time-step).
     */
    RunnableVector &getPreTimeStepRunnableVector()
    {
      return m_preRunnableVector;
    }

    /**
     * Adds reference to an object whose 'run' method is executed at
     * the end of the runOneStep method.
     */
    void addPostTimeStepRunnable(esys::lsm::Runnable &runnable);

    /**
     * Returns vector of Runnable objects (which are run post-time-step).
     */
    const RunnableVector &getPostTimeStepRunnableVector() const
    {
      return m_postRunnableVector;
    }

    /**
     * Returns vector of Runnable objects (which are run post-time-step).
     */
    RunnableVector &getPostTimeStepRunnableVector()
    {
      return m_postRunnableVector;
    }

    /**
     * Returns the dimensions relevant to calculating the minimum and maximum extents of
     * particles within a domain, taking into account periodic boundaries and the
     * dimensionality of the problem.
     */
    esys::lsm::IntVector getParticleDimensions()
    {
      return m_particle_dimensions;
    }

    /**
     * Provides the initial minimum and maximum extents of all the particles read in from a geometry file.
     *
     * @param initMinPt Minimum extent of particles inside domain.
     * @param initMaxPt Maximum extent of particles inside domain.
     */
    void getInitMinMaxPt(Vec3 &initMinPt, Vec3 &initMaxPt);

    /**
     * Defines the bounding box in which particles may roam.
     */
    void setSpatialDomain(const Vec3 &minBBoxPt, const Vec3 &maxBBoxPt);

    /**
     * Defines the bounding box in which particles may roam.
     */
    void setSpatialDomain(
      const Vec3 &minBBoxPt,
      const Vec3 &maxBBoxPt,
      const esys::lsm::IntVector &circDimVector
    );

    /**
     * Returns whether the setSpacialDomain method has been called with a
     * non-zero volume bounding box.
     */
    //bool haveSetSpatialDomain() const;

    /**
     * Gathers slave dimensions and coordinates.
     */
    void getSlaveSpatialDomains();


    /**
     * Template method for reading geometry from a specified file.
     * The template parameter is the type of particles created from
     * the particle data contained in the geometry file.
     *
     * @param fileName the name of the file containing geometry info.
     */
    template <class TmplParticle>
    void  readGeometry(const std::string &fileName);

    /**
     * Reads an initial particle configuration from a geometry file.
     *
     * @param fileName Name of file containing geometry information.
     */
    void readGeometryFile(const std::string &fileName);

    /**
     * Reads particle and interaction configuration from a check-point
     * summary file.
     *
     * @param fileName Name of file containing summary of check-point
     *                 data.
     */
    void loadCheckPointData(const std::string &checkPointFileName);

    /**
     * Adds enumeration of particles to the lattice.
     *
     * @param it Iterator object with next() and hasNext() methods.
     *           The TmplIterator::next() method is required to
     *           return an object which can accepted as a constructor
     *           argument for the CParticle class.
     */
    template <class TmplIterator, class TmplParticle>
    void addParticles(TmplIterator &it);

    /**
     * Adds enumeration of connections/bonds to the lattice.
     *
     * @param it Iterator object with next() and hasNext() methods.
     *           The TmplIterator::next() method is required to
     *           return an object which can accepted as a constructor
     *           argument for the SimpleConnectionData class.
     */
    template <class TmplIterator>
    void addConnections(TmplIterator &it);

    //--- function for mesh data exchange ---
    template <typename TmplVisitor>
    void visitMeshFaceReferences(const string &meshName);

    template <typename TmplVisitor>
    void visitMesh2dNodeReferences(const string &meshName, TmplVisitor &visitor);

    template <typename TmplVisitor>
    void visitMesh2dEdgeStress(const string &meshName, TmplVisitor &visitor);

    template <typename TmplVisitor>
    void visitTriMeshFaceForce(
      const string &meshName,
      TmplVisitor &visitor
    );

    typedef std::vector<int> IdVector;

    template <typename TmplVisitor, typename TmplParticle>
    void visitParticlesOfType(
      const IdVector &particleIdVector,
      TmplVisitor &visitor
    );

    template <typename TmplVisitor>
    void visitParticles(const IdVector &particleIdVector, TmplVisitor &visitor);

	// --- console related fucntions
    void setVerbosity(int);
	void initializeConsole(const string&, int);
	void setConsoleFilename(const string&);
	void setConsoleBuffered(unsigned int);
	
	// --- setting interaction parameters during a simulation ---
	void setInteractionParameter(const string&,const string&,double);

protected:
    /*!
     * Updates the minimum and maximum extents of the particles as each is read in.
     *
     * @param particle Information (e.g., position, radius) for a particle.
     *
     */
    template<typename TmplParticle>
    void particlesMinMax(const TmplParticle &particle);
};

#include "Parallel/LatticeMaster.hpp"

#endif

