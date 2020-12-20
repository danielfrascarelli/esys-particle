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

#include "Parallel/LatticeMaster.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

//--- IO includes --
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>

using std::ifstream;
using std::ofstream;
using std::fstream;
using std::ios;

// -- STL includes --
#include <map>
#include <algorithm>

using std::map;

// -- project includes --
#include "Model/TempPartStore.h"
#include "Parallel/mpi_tag_defs.h"
#include "Parallel/mpicmdbuf.h"
#include "Parallel/mpisgvbuf.h"
#include "Parallel/sublattice_cmd.h"
#include "Parallel/mpibarrier.h"
#include "Parallel/BroadCast_cmd.h"
#include "Parallel/BMesh2D_cmd.h"
#include "Model/Damping.h"
#include "Model/LocalDamping.h"
#include "Model/ElasticInteractionGroup.h"
#include "Model/Wall.h"
#include "Model/BWallInteractionGroup.h"
#include "Model/ViscWallIG.h"
#include "Model/SoftBWallInteractionGroup.h"
#include "Model/BodyForceGroup.h"
#include "Foundation/Timer.h"
#include "Foundation/PathSearcher.h"
#include "Foundation/StringUtil.h"
#include "Foundation/Functional.h"
#include "Fields/ParticleFieldMaster.h"
#include "Fields/ScalarParticleDistributionMaster.h"
#include "Fields/InteractionFieldMaster.h"
#include "Fields/VectorInteractionFieldMaster.h"
#include "Fields/VectorTriangleFieldMaster.h"
#include "Fields/ScalarTriangleFieldMaster.h"
#include "Fields/WallFieldMaster.h"
#include "Fields/TriggeredVectorParticleFieldMaster.h"
#include "Parallel/CheckPointInfo.h"
#include "Parallel/CheckPointLoader.h"
#include "Parallel/CheckPointController.h"
#include "Parallel/GeometryReader.h"
#include "Foundation/BoundingBox.h"
#include "Parallel/MeshReader.h"
#include "Parallel/Mesh2DReader.h"
#include "Model/MeshData2D.h"
#include "Parallel/MpiInfo.h"

#include "Parallel/MpiWrap.h"

#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include <boost/filesystem/path.hpp>
#include <cstdlib>

using namespace esys::lsm;

//std::string CLatticeMaster::s_workerExeName="ESySParticleBEWorker";

CLatticeMaster::CLatticeMaster()
  : m_timingFileName(),
    m_pTimers(new MpiWTimers()),
    m_pCheckPointController(new CheckPointController()),
    m_pSnapShotController(new CheckPointController()),
    m_processDims(3, 0),
    m_save_fields(),
    m_bbx_has_been_set(false),
    m_geometry_is_initialized(false),
    m_global_rank(0),
    m_global_size(0),
    m_max_ts(0),
    m_center_id(0),
    m_total_time(0),
    m_t(0),
    m_t_f(0), //fluid contents
    m_t_old(0), //fluid contents
    m_dt(0),
    m_dt_f(0), //fluid contents
    m_isInitialized(false),
    m_fluidinitiated(false), //fluid contents
    m_particle_type(),
    m_preRunnableVector(),
    m_postRunnableVector(),
    m_tml_global_comm(MPI_COMM_NULL),
    m_global_comm(MPI_COMM_NULL),
    m_dbl_NaN(std::numeric_limits<double>::quiet_NaN()),
    m_init_min_pt(m_dbl_NaN,m_dbl_NaN,m_dbl_NaN),
    m_init_max_pt(m_dbl_NaN,m_dbl_NaN,m_dbl_NaN),
    m_particle_dimensions(3)
{
  m_first_time=true;
}

CLatticeMaster::~CLatticeMaster()
{
  console.Debug() << "CLatticeMaster::~CLatticeMaster: enter\n";
  delete m_pTimers;
  delete m_pCheckPointController;
  delete m_pSnapShotController;
  for (vector<AFieldMaster*>::iterator it = m_save_fields.begin(); it != m_save_fields.end(); it++)
  {
    delete (*it);
  }
  m_save_fields.clear();

  // cleanup MPI groups and communicators
  // MPI_Group_free(&m_mpi_global_group);
  MPI_Group_free(&m_mpi_local_group);

  //MPI_Comm_free(&m_local_comm);

  console.Debug() << "CLatticeMaster::~CLatticeMaster: exit\n";
}

/****fluid contents: begin****/
/*!
  add fluid to the sublattice

  \param Bw Bulk modulus of water
  \param Bp Bulk modulus of particle
  \param Mu Viscosity of water
  \param alpha Adjusting factor between two time steps
*/
void CLatticeMaster::addFluidInteraction(double cellside,double Bw, double Bp, double Mu, double alpha, double flowrate, double pressure, Vec3 inflow, Vec3 outflow, double dt_f)
{
  console.Debug() << "CLatticeMaster::addFluidInteraction\n";

  BroadcastCommand cmd(getGlobalRankAndComm(), CMD_ADDFLUID);

  cmd.append(cellside);
  cmd.append(Bw);
  cmd.append(Bp);
  cmd.append(Mu);
  cmd.append(alpha);
  cmd.append(flowrate);
  cmd.append(pressure);
  cmd.append(inflow);
  cmd.append(outflow);
  //send command to slave
  cmd.broadcast();
  m_dt_f=dt_f;
  m_fluidinitiated=true;
  console.Debug() << "end CLatticeMaster::addFluidInteraction() \n";
}

/*!
  add fluid to the sublattice

  \param Bw Bulk modulus of water
  \param Bp Bulk modulus of particle
  \param Mu Viscosity of water
  \param alpha Adjusting factor between two time steps
*/
void CLatticeMaster::addFluidInteractionVec3(Vec3 cellside,double Bw, double Bp, double Mu, double alpha, double flowrate, double pressure, Vec3 inflow, Vec3 outflow, double dt_f)
{
  console.Debug() << "CLatticeMaster::addFluidInteractionVec3\n";

  BroadcastCommand cmd(getGlobalRankAndComm(), CMD_ADDFLUID_VEC3);

  cmd.append(cellside);
  cmd.append(Bw);
  cmd.append(Bp);
  cmd.append(Mu);
  cmd.append(alpha);
  cmd.append(flowrate);
  cmd.append(pressure);
  cmd.append(inflow);
  cmd.append(outflow);
  //send command to slave
  cmd.broadcast();
  m_dt_f=dt_f;
  m_fluidinitiated=true;
  console.Debug() << "end CLatticeMaster::addFluidInteractionVec3() \n";
}


/*!
  add a scalar fluid field to the list of fields to be saved

  \param filename the name of the file the field is saved into
  \param fieldname the name of the field
  \param savetype output file format
  \param t_0 first timestep to be saved
  \param t_end last timestep to be saved
  \param dt timesteps between saves
*/
void CLatticeMaster::addScalarFluidSaveField(
  const string& filename,
  const string& fieldname,
  const string& savetype,
  int t_0,
  int t_end,
  int dt
)
{
  console.Debug()
    << "CLatticeMaster::addScalarFluidSaveField("
    << filename
    << ","
    << fieldname
    << ","
    << t_0
    << ","
    << t_end
    << ","
    << dt
    << ")\n";

  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CMPIBarrier barrier(m_global_comm);

  //send cmd to slave
  cmd_buffer.broadcast(CMD_ADD_SFF);

  AFieldMaster* new_fm =
    new ScalarFluidFieldMaster(
      &m_tml_global_comm,
      fieldname,
      filename,
      savetype,
      t_0,
      t_end,
      dt
    );
  m_save_fields.push_back(new_fm);

  barrier.wait("CLatticeMaster::addScalarFluidSaveField");

  console.Debug() << "end CLatticeMaster::addScalarFluidSaveField()\n";
}

/*!
  add a vector fluid field to the list of fields to be saved

  \param filename the name of the file the field is saved into
  \param fieldname the name of the field
  \param savetype output file format
  \param t_0 first timestep to be saved
  \param t_end last timestep to be saved
  \param dt timesteps between saves
*/
void CLatticeMaster::addVectorFluidSaveField(
  const string& filename,
  const string& fieldname,
  const string& savetype,
  int t_0,
  int t_end,
  int dt
)
{
  console.Debug()
    << "CLatticeMaster::addVectorFluidSaveField("
    << filename
    << ","
    << fieldname
    << ","
    << t_0
    << ","
    << t_end
    << ","
    << dt
    << ")\n";

  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CMPIBarrier barrier(m_global_comm);

  //send cmd to slave
  cmd_buffer.broadcast(CMD_ADD_VFF);

  AFieldMaster* new_fm =
    new VectorFluidFieldMaster(
      &m_tml_global_comm,
      fieldname,
      filename,
      savetype,
      t_0,
      t_end,
      dt
    );
  m_save_fields.push_back(new_fm);

  barrier.wait("CLatticeMaster::addVectorFluidSaveField");

  console.Debug() << "end CLatticeMaster::addVectorFluidSaveField()\n";
}

/*!
  add a scalar fluid interaction field to the list of fields to be saved

  \param filename the name of the file the field is saved into
  \param fieldname the name of the field
  \param savetype output file format
  \param t_0 first timestep to be saved
  \param t_end last timestep to be saved
  \param dt timesteps between saves
*/
void CLatticeMaster::addScalarFluidInteractionSaveField(
  const string& filename,
  const string& fieldname,
  const string& savetype,
  int t_0,
  int t_end,
  int dt
)
{
  console.Debug()
     << "CLatticeMaster::addScalarFluidInteractionSaveField("
     << filename
     << ","
     << fieldname
     << ","
     << t_0
     << ","
     << t_end
     << ","
     << dt
     << ")\n";

  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CMPIBarrier barrier(m_global_comm);

  //send cmd to slave
  cmd_buffer.broadcast(CMD_ADD_SFIF);

  AFieldMaster* new_fm =
    new ScalarFluidInteractionFieldMaster(
      &m_tml_global_comm,
      fieldname,
      filename,
      savetype,
      t_0,
      t_end,
      dt
    );
  m_save_fields.push_back(new_fm);

  barrier.wait("CLatticeMaster::addScalarFluiedInteractionSaveField");

  console.Debug() << "end CLatticeMaster::addScalarFluidInteractionSaveField()\n";
}

/*!
  add a vector fluid interaction field to the list of fields to be saved

  \param filename the name of the file the field is saved into
  \param fieldname the name of the field
  \param savetype output file format
  \param t_0 first timestep to be saved
  \param t_end last timestep to be saved
  \param dt timesteps between saves
*/
void CLatticeMaster::addVectorFluidInteractionSaveField(
  const string& filename,
  const string& fieldname,
  const string& savetype,
  int t_0,
  int t_end,
  int dt
)
{
  console.Debug()
     << "CLatticeMaster::addVectorFluidInteractionSaveField("
     << filename
     << ","
     << fieldname
     << ","
     << t_0
     << ","
     << t_end
     << ","
     << dt
     << ")\n";

  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CMPIBarrier barrier(m_global_comm);

  //send cmd to slave
  cmd_buffer.broadcast(CMD_ADD_VFIF);

  AFieldMaster* new_fm =
    new VectorFluidInteractionFieldMaster(
      &m_tml_global_comm,
      fieldname,
      filename,
      savetype,
      t_0,
      t_end,
      dt
    );
  m_save_fields.push_back(new_fm);

  barrier.wait("CLatticeMaster::addVectorFluiedInteractionSaveField");

  console.Debug() << "end CLatticeMaster::addVectorFluidInteractionSaveField()\n";
}

/****fluid contents: end****/


void CLatticeMaster::init()
{
  MPI_Comm_size(MPI_COMM_WORLD, &m_global_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &m_global_rank);

  console.Debug() << "Master has rank " << m_global_rank << "\n" ;
  console.Debug() << "Global size is  " << m_global_size << "\n" ;
}

void CLatticeMaster::setupWorkers(int ns)
{
  // setup communicators and check size
  m_global_comm=MPI_COMM_WORLD;
  MPI_Comm_size(m_global_comm, &m_global_size);
  m_tml_global_comm.setComm(m_global_comm);

  // -- local communicator = global comm - proc(0)
  // get global MPI_Group
  MPI_Group mpi_global_group;
  MPI_Comm_group(MPI_COMM_WORLD,&mpi_global_group);
  // subtract id 0 from global group
  int id0=0;
  MPI_Group_excl(mpi_global_group,1,&id0,&m_mpi_local_group);
  //  create communicator
  MPI_Comm_create(MPI_COMM_WORLD,m_mpi_local_group,&m_local_comm);

  // free global group
  MPI_Group_free(&mpi_global_group);

  m_tml_global_comm.barrier();

  if(m_global_size!=ns+1){
    cerr << "wrong number of processes !! aborting" << endl << flush;
    exit(255);
    // should send abort command to workers
  }
  initializeConsole("console.out",1000); // 1st barrier -> sync time
}

int CLatticeMaster::getNumWorkerProcesses() const
{
  return m_global_size - 1;
}

void CLatticeMaster::setNumSteps(int s)
{
  m_max_ts=s;
  m_pCheckPointController->setNumTimeSteps(s);
  m_pSnapShotController->setNumTimeSteps(s);
}

void CLatticeMaster::saveTimingDataToFile(const std::string &fileNamePrefix)
{
  std::stringstream sStream;
  sStream << fileNamePrefix << m_global_rank << ".csv";
  setTimingFileName(sStream.str());

  CMPILCmdBuffer cmdBuffer(m_global_comm, MpiInfo(m_global_comm).rank());
  CMPIBarrier barrier(m_global_comm);
  cmdBuffer.broadcast(CMD_PERFORMTIMING);

  CVarMPIBuffer buffer(m_global_comm);
  // send timing filename prefix
  buffer.append(fileNamePrefix.c_str());
  buffer.broadcast(0);
  buffer.clear();
  m_tml_global_comm.barrier("performTiming");
}

void CLatticeMaster::setTimingFileName(const std::string &fileName)
{
  m_timingFileName = fileName;
}

const std::string &CLatticeMaster::getTimingFileName() const
{
  return m_timingFileName;
}

void CLatticeMaster::do2dCalculations(bool do2d)
{
  // Set the 2-D or 3-D dimension information for checkpoint files.
  m_pCheckPointController->set_is2d(do2d);
  m_pSnapShotController->set_is2d(do2d);

  CMPILCmdBuffer cmdBuffer(m_global_comm, MpiInfo(m_global_comm).rank());
  CMPIBarrier barrier(m_global_comm);
  cmdBuffer.broadcast(CMD_DO2DCALCULATIONS);

  CVarMPIBuffer buffer(m_global_comm);
  buffer.append(static_cast<int>(do2d));
  buffer.broadcast(0);

  m_tml_global_comm.barrier("do2dCalculations cmd");
}

/*!
  Make a lattice from particles of the given type. Does not set up neighbor tables.

  \param ptype the type of particle
  \param nrange search range
  \param alpha pair search cutoff
*/
void CLatticeMaster::makeLattice(const char* ptype, double nrange, double alpha)
{
  console.XDebug()
     << "CLatticeMaster::makeLattice( ptype=" << ptype
     << " , nrange=" << nrange << ", alpha=" << alpha << " )\n";

  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer buffer(m_global_comm);

  m_particle_type = string(ptype);

  //send command to slave
  cmd_buffer.broadcast(CMD_MAKELATTICE);

  // send lattice parameters

  CLatticeParam clp(ptype, nrange, alpha, getProcessDims());
  clp.packInto(&buffer);
  buffer.broadcast(0);
  buffer.clear();

  m_tml_global_comm.barrier("CLatticeMaster::makeLattice");
  //m_tml_global_comm.barrier();

  console.XDebug() << "end CLatticeMaster::makeLattice\n";
}

/*!
  make lattice and set time step size

  \param ptype the type of particle
  \param nrange search range
  \param alpha pair search cutoff
  \param dt time step
*/
void CLatticeMaster::makeLattice(const char* ptype, double nrange, double alpha, double dt)
{
  makeLattice(ptype, nrange, alpha);
  setTimeStepSize(dt);
}

/*!
  set the time step size

  \param dt time step
*/
void CLatticeMaster::setTimeStepSize(double dt)
{
  m_pCheckPointController->setTimeStepSize(dt);
  m_pSnapShotController->setTimeStepSize(dt);
  m_dt = dt;

  BroadcastCommand cmd(getGlobalRankAndComm(), CMD_SETTIMESTEPSIZE);
  cmd.append(m_dt);
  cmd.broadcast();
}

void CLatticeMaster::setProcessDims(const CLatticeParam::ProcessDims &dims)
{
  m_processDims = dims;
}

const CLatticeParam::ProcessDims &CLatticeMaster::getProcessDims() const
{
  return m_processDims;
}

/*!
 * Provides the initial minimum and maximum extents of all the particles read in from a geometry file.
 *
 * \param initMinPt Minimum extent of particles inside domain.
 * \param initMaxPt Maximum extent of particles inside domain.
 */
void CLatticeMaster::getInitMinMaxPt(Vec3 &initMinPt, Vec3 &initMaxPt)
{
  initMinPt = m_init_min_pt;
  initMaxPt = m_init_max_pt;
}

/*!
   Define model bounding box and setup neighbor tables in Workers. Non-circular version.

   \param minBBoxPt minimum point of the bounding box
   \param maxBBoxPt maximum point of the bounding box
*/
void CLatticeMaster::setSpatialDomain(const Vec3 &minBBoxPt, const Vec3 &maxBBoxPt)
{
  console.Debug() << "CLatticeMaster::setSpatialDomain: enter\n";
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);

  // check if domain has already been set
  if(m_geometry_is_initialized){
    console.Warning() << "Spatial Domain has already been set and can not be changed - this call does nothing.\n";
  } else {
    // set geometry info according to parameters
    // no need to set circular bc - defaults to (0,0,0) anyway
    m_geo_info.setBBox(minBBoxPt,maxBBoxPt);
    m_bbx_has_been_set=true;
    // none of the boundaries circular -> simple initLattice
    cmd_buffer.broadcast(CMD_INITLATTICE);
    esys::lsm::Vec3Vector corners;
    corners.push_back(minBBoxPt);
    corners.push_back(maxBBoxPt);
    // send data to workers
    m_tml_global_comm.broadcast_cont_packed(corners);
    getSlaveSpatialDomains();
    m_tml_global_comm.barrier(std::string("cmd=") + StringUtil::toString(CMD_INITLATTICE));
    m_geometry_is_initialized=true;
  }
  console.Debug() << "CLatticeMaster::setSpatialDomain: exit\n";
}

/*!
   Define model bounding box and setup neighbor tables in Workers. Circular boundary version.

   \param minBBoxPt minimum point of the bounding box
   \param maxBBoxPt maximum point of the bounding box
   \param circDimVector a vector of ints containing the circular boundary condition flags for all dimensions
*/
void CLatticeMaster::setSpatialDomain(
  const Vec3 &minBBoxPt,
  const Vec3 &maxBBoxPt,
  const esys::lsm::IntVector &circDimVector
)
{
  console.Debug() << "CLatticeMaster::setSpatialDomain circular: enter\n";
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  esys::lsm::BoolVector periodicDimensions;

  // check if domain has already been set
  if(m_geometry_is_initialized){
    console.Warning() << "Spatial Domain has already been set and can not be changed - this call does nothing.\n";
  } else {
    // check if we actually need the circular part or if the periodic dimensions are (false,false,false)
    if((circDimVector[0]!=0) || (circDimVector[1]!=0) || (circDimVector[2]!=0)) { // we have at least one periodic dimension
        // set geometry info according to parameters
        m_geo_info.setBBox(minBBoxPt,maxBBoxPt);
        // need to set circular bc
        for (int i = 0; i <= 2; i++) periodicDimensions.push_back(bool(circDimVector[i]));
        m_geo_info.setPeriodicDimensions(periodicDimensions);
        m_bbx_has_been_set=true;

        // at least one boundary circular -> initLatticeCirc
        console.XDebug()<< "Circular Lattice Init\n";
        cmd_buffer.broadcast(CMD_INITLATTICECIRC);

        esys::lsm::Vec3Vector corners;
        corners.push_back(minBBoxPt);
        corners.push_back(maxBBoxPt);

        console.XDebug()<< "Broadcasting bounding-box corners\n";
        m_tml_global_comm.broadcast_cont_packed(corners);

        // send circ. markers
        m_tml_global_comm.broadcast_cont(circDimVector);

        getSlaveSpatialDomains();
        m_tml_global_comm.barrier("CMD_INITLATTICECIRC");
        console.Debug() << "CLatticeMaster::setSpatialDomain circular: exit\n";
        m_geometry_is_initialized=true;
        m_bbx_has_been_set=true;
    } else { // all dimensions non-periodic -> call non-circular version of setSpatialDomain
        setSpatialDomain(minBBoxPt,maxBBoxPt);
    }
   }
}



void CLatticeMaster::getSlaveSpatialDomains()
{
  console.XDebug()
    <<"CLatticeMaster::getSlaveSpatialDomains: enter\n";

  // get slave dimensions and coordinates
  multimap<int,int> slavedims;
  multimap<int,int> slavecoords;
  console.XDebug()
    <<"CLatticeMaster::getSlaveSpatialDomains: gathering dims\n";

  m_tml_global_comm.gather(slavedims);

  console.XDebug()
    <<"CLatticeMaster::getSlaveSpatialDomains: gathering coords\n";
  m_tml_global_comm.gather(slavecoords);

  console.XDebug()
    <<"CLatticeMaster::getSlaveSpatialDomains: received dims & coords\n";

  // extract slave dimensions from received data
  int dim_x,dim_y,dim_z; // dimensions in x,y,z direction
  multimap<int,int>::iterator iter=slavedims.find(1);
  dim_x=iter->second;
  iter++;
  dim_y=iter->second;
  iter++;
  dim_z=iter->second;
  // setup coord->rank mapping for slaves
  int nslaves=dim_x*dim_y*dim_z; // number of slaves
  for(int i=0;i<nslaves;i++){
    // extract slave coordinates from received data
    //int cx,cy,cz; // x-,y-,z- coord.

    multimap<int,int>::iterator coord_iter=slavecoords.find(i+1); // slaves are 1..n
    //cx=coord_iter->second;
    coord_iter++;
    //cy=coord_iter->second;
    coord_iter++;
    //cz=coord_iter->second;
  }
  console.XDebug()
    <<"CLatticeMaster::getSlaveSpatialDomains: exit\n";

}

/*!
  read geometry file
  calls appropriate reader function depending on particle type

  \param fileName the name of the geometry file
*/
void  CLatticeMaster::readGeometryFile(const string& fileName)
{
  console.Debug()<<"CLatticeMaster::readGeometryFile( " << fileName << ")\n";

  if (m_particle_type == "Basic"){
    readGeometry<CParticle>(fileName);
  } else if (m_particle_type == "Rot") {
    readGeometry<CRotParticle>(fileName);
  } else if (m_particle_type == "RotVi") {
    readGeometry<CRotParticleVi>(fileName);
  } else if (m_particle_type == "RotTherm") {
    readGeometry<CRotThermParticle>(fileName);
  } else {
    throw std::runtime_error(
      std::string("Unknown particle type: ")
      +
      m_particle_type
    );
  }
}

/*!
  Load data from a save checkpoint in order to restart the simulation.
  Control flow:
  CLatticeMaster::loadCheckPointData() (i.e. here)
  ->  CheckPointController::issueCheckPointLoadingCmd()
  ->  [on worker] CheckPointer::loadCheckPoint()
  ->  SubLattice::loadCheckPointData()

  \param checkPointFileName the name of the base file (*_0.txt) of the checkpoint
*/
void CLatticeMaster::loadCheckPointData(const string &checkPointFileName)
{
  CheckPointInfo cpInfo;

  std::ifstream iStream(checkPointFileName.c_str());
  if (!iStream) {
    throw std::runtime_error(
      std::string("Could not open file ") + checkPointFileName
    );
  }

  // read geometry info from checkpoint file and initialise geometry
  cpInfo.read(iStream);
  GeometryInfo geoInfo = cpInfo.getGeometryInfo();

  // check if geometry info has been set previously
  if(m_bbx_has_been_set){
    // check if old & new geometry are identical (the "compatible" check of readGeometry is not sufficient)
    if(!m_geo_info.isIdenticalGeometry(geoInfo)){
    std::cerr << "Geometry info read from file is different from previously set geometry (bounding box, circular boundaries) - Model may not run properly!" << std::endl;
    }
  } else {
    m_geo_info=geoInfo; // set geometry meta-info

    // set model spatial domain -> sets geometry info & initializes neighbor table in workers
    if(m_geo_info.hasAnyPeriodicDimensions()){
    setSpatialDomain(m_geo_info.getBBoxCorners()[0],m_geo_info.getBBoxCorners()[1],m_geo_info.getPeriodicDimensions());
    } else {
    setSpatialDomain(m_geo_info.getBBoxCorners()[0],m_geo_info.getBBoxCorners()[1]);
    }
    m_bbx_has_been_set=true; // flag that spatial domain has been set
  }

  m_t=cpInfo.getTimeStep();
  std::cout << "set time step to " << m_t << std::endl;

  // read data
  m_pCheckPointController->issueCheckPointLoadingCmd(checkPointFileName);
  searchNeighbors(true);
  updateInteractions();
}

/*!
  add a wall to the sublattice

  \param wname the name of the wall
  \param ipos initial position of the wall
  \param inorm initial normal of the wall
*/
void CLatticeMaster::addWall(const string& wname, const Vec3& ipos,const Vec3& inorm)
{
  console.Debug() << "CLatticeMaster::addWall("<< wname << " , " << ipos << " , " << inorm << ")\n" ;

  BroadcastCommand cmd(getGlobalRankAndComm(), CMD_ADDWALL);

  cmd.append(wname.c_str());
  cmd.append(ipos);
  cmd.append(inorm);
  //send command to slave
  cmd.broadcast();

  console.Debug() << "end CLatticeMaster::initWall() \n";
}



/*!
  add an elastic wall interaction to model

  \param param the interaction parameters (name of the wall, name of the interaction, spring constant)
*/
void CLatticeMaster::addWallIG(const CEWallIGP& param)
{
  console.Debug() << "CLatticeMaster::addWallIG( -elastic- " << param.getWallName() << " , " << param.getName() << " , " << param.getSpringConst() << ")\n" ;

  BroadcastCommand cmd(getGlobalRankAndComm(), CMD_ADDEWALLIG);

  //send command to slave
  cmd.packInto(param);
  //send command to slave
  cmd.broadcast();

  console.Debug() << "end CLatticeMaster::addWallIG() \n";
}

/*!
  add a bonded wall interaction to model

  \param param the interaction parameters (name of the wall, name of the interaction, spring constant)
*/
void CLatticeMaster::addWallIG(const CBWallIGP& param)
{
  console.Debug() << "CLatticeMaster::addWallIG( -bonded- "
          << param.getWallName() << " , "
          << param.getName() << " , "
          << param.getSpringConst() << " , "
          << param.getTag() << ")\n" ;

  BroadcastCommand cmd(getGlobalRankAndComm(), CMD_ADDBWALLIG);

  //send command to slave
  cmd.packInto(param);
  //send command to slave
  cmd.broadcast();

  console.Debug() << "end CLatticeMaster::addWallIG() \n";
}

/*!
  add a viscous wall interaction to model

  \param param the interaction parameters (name of the wall, name of the interaction, spring constant)
*/
void CLatticeMaster::addWallIG(const CVWallIGP& param)
{
  console.Debug() << "CLatticeMaster::addWallIG( -viscous- "
          << param.getWallName() << " , "
          << param.getName() << " , "
          << param.getSpringConst() << " , "
          << param.getNu() << " , "
          << param.getTag() << ")\n" ;

  BroadcastCommand cmd(getGlobalRankAndComm(), CMD_ADDBWALLIG);

  //send command to slave
  cmd.packInto(param);
  //send command to slave
  cmd.broadcast();

  console.Debug() << "end CLatticeMaster::addWallIG() \n";
}
/*!
  add a bonded wall interaction with direction-dependent spring constant to model

  \param param the interaction parameters (name of the wall, name of the interaction, spring constant)
*/
void CLatticeMaster::addWallIG(const CSoftBWallIGP& param)
{
  console.Debug() << "CLatticeMaster::addWallIG( -directional bonded- "
          << param.getWallName() << " , "
          << param.getName() << " , "
          << param.getNormalK() << " , "
          << param.getShearK() << " , "
          << param.getTag() << " , "
          << param.getScaling() << ")\n" ;

  BroadcastCommand cmd(getGlobalRankAndComm(), CMD_ADDBBWALLIG);

  //send command to slave
  cmd.packInto(param);
  //send command to slave
  cmd.broadcast();

  console.Debug() << "end CLatticeMaster::addWallIG() \n";
}

/*!
  add a tagged elastic wall interaction to the model

  \param param the interaction parameters (name of the wall, name of the interaction, spring constant)
  \param tag the tag of the particles the wall is interacting with
  \param mask the mask determining which bits of the tag are significant
*/
void CLatticeMaster::addTaggedWallIG(const CEWallIGP& param,int tag,int mask)
{
  console.Debug() << "CLatticeMaster::addTaggedWallIG( -elastic- " << param.getWallName() << " , " << param.getName() << " , " << param.getSpringConst() << ")\n" ;

  BroadcastCommand cmd(getGlobalRankAndComm(), CMD_ADDTAGGEDEWALLIG);

  //pack data into command buffer
  cmd.packInto(param);
  cmd.append(tag);
  cmd.append(mask);
  //send command to slave
  cmd.broadcast();

  console.Debug() << "end CLatticeMaster::addTaggedWallIG() \n";
}

/*!
  Get position of a wall. Returns (0,0,0) for a non-existing wall.

  \param WallName the name of the wall
*/
Vec3 CLatticeMaster::getWallPosn(const std::string& WallName)
{
  console.Debug() << "CLatticeMaster::getWallPosn: enter\n";
  Vec3 wpos=Vec3(0.0,0.0,0.0);

  // setup command
  BroadcastCommand cmd(getGlobalRankAndComm(), CMD_GETWALLPOS);

  // pack name into command buffer
  cmd.append(WallName.c_str());

  // broadcast command buffer
  cmd.broadcast();

  // collect data
  multimap<int,Vec3> posmap;
  m_tml_global_comm.gather(posmap);

  wpos=(posmap.find(1))->second;

  console.Debug() << "CLatticeMaster::getWallPosn: exit\n";
  return wpos;
}

/*!
  Get force acting on a wall. Returns (0,0,0) for a non-existing wall.

  \param WallName the name of the wall
*/
Vec3 CLatticeMaster::getWallForce(const std::string& WallName)
{
  console.Debug() << "CLatticeMaster::getWallForce: enter\n";
  Vec3 wforce=Vec3(0.0,0.0,0.0);

  // setup command
  BroadcastCommand cmd(getGlobalRankAndComm(), CMD_GETWALLFORCE);

  // pack name into command buffer
  cmd.append(WallName.c_str());

  // broadcast command buffer
  cmd.broadcast();

  // collect data
  multimap<int,Vec3> forcemap;
  m_tml_global_comm.gather(forcemap);

  for(multimap<int,Vec3>::iterator iter=forcemap.begin();
      iter!=forcemap.end();
      iter++){
    wforce=wforce+iter->second;
  }

  console.Debug() << "CLatticeMaster::getWallForce: exit\n";
  return wforce;
}

/*!
  add a sphere body to the sublattice

  \param sname the name of the sphere
  \param ipos initial position of the sphere
  \param radius radius of the sphere
*/
void CLatticeMaster::addSphereBody(const string& sname, const Vec3& ipos,const double& radius)
{
  console.Debug() << "CLatticeMaster::addSphereBody("<< sname << " , " << ipos << " , " << radius << ")\n" ;

  BroadcastCommand cmd(getGlobalRankAndComm(), CMD_ADDSPHEREBODY);

  cmd.append(sname.c_str());
  cmd.append(ipos);
  cmd.append(radius);
  //send command to slave
  cmd.broadcast();

  console.Debug() << "end CLatticeMaster::addSphereBody() \n";
}

/*!
  add an elastic sphere body interaction to model

  \param param the interaction parameters (name of the sphere, name of the interaction, spring constant)
*/
void CLatticeMaster::addSphereBodyIG(const CESphereBodyIGP& param)
{
  console.Debug() << "CLatticeMaster::addSphereBodyIG( -elastic- " << param.getSphereBodyName() << " , " << param.getName() << " , " << param.getSpringConst() << ")\n" ;

  BroadcastCommand cmd(getGlobalRankAndComm(), CMD_ADDESPHEREBODYIG);

  //send command to slave
  cmd.packInto(param);
  //send command to slave
  cmd.broadcast();

  console.Debug() << "end CLatticeMaster::addSphereBodyIG() \n";
}

/*!
  Get position of a sphere body. Returns (0,0,0) for a non-existing sphere body.

  \param SphereName the name of the wall
*/
Vec3 CLatticeMaster::getSphereBodyPosn(const std::string& SphereName)
{
  console.Debug() << "CLatticeMaster::getSphereBodyPosn: enter\n";
  Vec3 wpos=Vec3(0.0,0.0,0.0);

  // setup command
  BroadcastCommand cmd(getGlobalRankAndComm(), CMD_GETSPHEREBODYPOS);

  // pack name into command buffer
  cmd.append(SphereName.c_str());

  // broadcast command buffer
  cmd.broadcast();

  // collect data
  multimap<int,Vec3> posmap;
  m_tml_global_comm.gather(posmap);

  wpos=(posmap.find(1))->second;

  console.Debug() << "CLatticeMaster::getSphereBodyPosn: exit\n";
  return wpos;
}

/*!
  Get force acting on a sphere body. Returns (0,0,0) for a non-existing sphere body.

  \param SphereName the name of the sphere body
*/
Vec3 CLatticeMaster::getSphereBodyForce(const std::string& SphereName)
{
  console.Debug() << "CLatticeMaster::getSphereBodyForce: enter\n";
  Vec3 wforce=Vec3(0.0,0.0,0.0);

  // setup command
  BroadcastCommand cmd(getGlobalRankAndComm(), CMD_GETSPHEREBODYFORCE);

  // pack name into command buffer
  cmd.append(SphereName.c_str());

  // broadcast command buffer
  cmd.broadcast();

  // collect data
  multimap<int,Vec3> forcemap;
  m_tml_global_comm.gather(forcemap);

  for(multimap<int,Vec3>::iterator iter=forcemap.begin();
      iter!=forcemap.end();
      iter++){
    wforce=wforce+iter->second;
  }

  console.Debug() << "CLatticeMaster::getSphereBodyForce: exit\n";
  return wforce;
}

/*!
  Change radius of all particles with a given tag by a given amount 

  \param tag the tag of the particles to be moved
  \param deltaR the amount of radius change
*/
void CLatticeMaster::changeRadiusBy(int particleTag, const double& deltaR)
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer buffer(m_global_comm);
  CMPIBarrier barrier(m_global_comm);

  cmd_buffer.broadcast(CMD_CHGRADIUS);
  buffer.append(particleTag);
  buffer.append(deltaR) ;
  buffer.broadcast(m_global_rank);
  barrier.wait("changeRadiusBy");
}

/*!
  Move all particles with a given tag to a given displacement relative to their original position.

  \param tag the tag of the particles to be moved
  \param d the displacement
*/
void CLatticeMaster::moveParticleTo(int tag, const Vec3& d)
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer buffer(m_global_comm);
  CMPIBarrier barrier(m_global_comm);

  cmd_buffer.broadcast(CMD_PMOVE);
  buffer.append(tag);
  buffer.append(d) ;
  buffer.broadcast(m_global_rank);
  barrier.wait("moveParticleTo");
}

/*!
  Move all particles with a given tag by a given displacement relative
  to their current position.

  \param tag the tag of the particles to be moved
  \param d the displacement
*/
void CLatticeMaster::moveTaggedParticlesBy(int tag, const Vec3& d)
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer buffer(m_global_comm);
  CMPIBarrier barrier(m_global_comm);

  cmd_buffer.broadcast(CMD_PMOVETAGGEDBY);
  buffer.append(tag);
  buffer.append(d) ;
  buffer.broadcast(m_global_rank);
  barrier.wait("moveTaggedParticlesBy");
  // force ntable rebuild
  searchNeighbors(true);
}

void CLatticeMaster::moveSingleParticleTo(int particleId, const Vec3& posn)
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer buffer(m_global_comm);
  CMPIBarrier barrier(m_global_comm);

  cmd_buffer.broadcast(CMD_IDPARTICLEMOVE);
  buffer.append(particleId);
  buffer.append(posn) ;
  buffer.broadcast(m_global_rank);
  barrier.wait("moveSingleParticleTo");
}

/*!
  Move a node in a TriMeshIG by a given amount.

  \param tm_name the name of the TriMeshIG
  \param id the node id
  \param d the displacement
*/
void CLatticeMaster::moveSingleNodeBy(const string& name,int id,const Vec3& d)
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer buffer(m_global_comm);
  CMPIBarrier barrier(m_global_comm);

  cmd_buffer.broadcast(CMD_MOVENODE);
  buffer.append(name.c_str());
  buffer.append(id) ;
  buffer.append(d);
  buffer.broadcast(m_global_rank);
  barrier.wait("moveSingleNodeBy");
}

/*!
  Move all nodes with a given tag in a TriMeshIG by a given amount.

  \param tm_name the name of the TriMeshIG
  \param tag the tag
  \param d the displacement
*/
void CLatticeMaster::moveTaggedNodesBy(const string& name,int tag,const Vec3& d)
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer buffer(m_global_comm);
  CMPIBarrier barrier(m_global_comm);

  cmd_buffer.broadcast(CMD_MOVETAGGEDNODES);
  buffer.append(name.c_str());
  buffer.append(tag) ;
  buffer.append(d);
  buffer.broadcast(m_global_rank);
  barrier.wait("moveTaggedNodesBy");
}

/*!
  Move a whole mesh by a given amount.

  \param meshName the name of the mesh to be moved
  \param translation the vector by which the mesh is translated
*/
void CLatticeMaster::translateMeshBy(
  const string &meshName,
  const Vec3& translation
)
{
  console.Debug() << "CLatticeMaster::translateMeshBy(): enter \n";
  BroadcastCommand cmd(getGlobalRankAndComm(), CMD_TRANSLATEMESHBY);

  cmd.append(meshName.c_str());
  cmd.append(translation);
  //send command to slave
  cmd.broadcast();
  console.Debug() << "CLatticeMaster::translateMeshBy(): exit \n";
}

/*!
  Rigidly rotate a whole mesh around a given axis.

  \param meshName the name of the mesh to be moved
  \param origin a point on the axis
  \param axis the orientation of the axis around which the mesh is rotated
  \param angle the rotation angle in radians
*/
void CLatticeMaster::rotateMeshBy(
    const string &meshName,
    const Vec3& origin,
    const Vec3& axis,
    const double angle
)
{
    console.Debug() << "CLatticeMaster::rotateMeshBy(): enter \n";
    BroadcastCommand cmd(getGlobalRankAndComm(), CMD_ROTATEMESHBY);

    cmd.append(meshName.c_str());
    cmd.append(origin);
    cmd.append(axis.unit()); // force normalised axis vector
    cmd.append(angle);
    //send command to worker
    cmd.broadcast();
    console.Debug() << "CLatticeMaster::rotateMeshBy(): exit \n";
    }

/*!
	Apply pressure to a mesh interaction. EXPERIMENTAL!!
	 
	\param meshName the name of the mesh to be moved
	\param pressure  
*/
void CLatticeMaster::applyPressureToMeshInteraction(const std::string& meshName,double pressure)
{
	console.Debug() << "CLatticeMaster::applyPressureToMeshInteraction (" << meshName << " , " << pressure << " )\n";
	BroadcastCommand cmd(getGlobalRankAndComm(), CMD_APPLYPRESSURRETOMESH);

	cmd.append(meshName.c_str());
	cmd.append(pressure);
	//send command to worker
	cmd.broadcast();
	
	console.Debug() << "CLatticeMaster::applyPressureToMeshInteraction() : exit \n";
}


/*!
  Add a given tag to the particle closest to a given position. Only the bits
  set in the mask will be influenced, i.e. new_tag=(old_tag & !mask) | (tag & mask).

  \param tag the tag
  \param mask the mask
  \param pos the position
*/
void CLatticeMaster::tagParticleNearestTo(int tag,int mask,const Vec3& pos)
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer buffer(m_global_comm);
  CMPIBarrier barrier(m_global_comm);

  cmd_buffer.broadcast(CMD_PTAG);
  buffer.append(tag);
  buffer.append(mask) ;
  buffer.append(pos);
  buffer.broadcast(m_global_rank);
  barrier.wait("tagParticleNearestTo");
}

/*!
  Add a given tag to the particle closest to a given position. Only the bits
  set in the mask will be influenced, i.e. new_tag=(old_tag & !mask) | (tag &
mask).

  \param tag the tag
  \param mask the mask
  \param pos the position
*/
int CLatticeMaster::findParticleNearestTo(const Vec3& pos)
{
  console.Debug() << "CLatticeMaster::findParticleNearestTo: enter\n";
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer buffer(m_global_comm);
  CMPIBarrier barrier(m_global_comm);

  cmd_buffer.broadcast(CMD_FINDNEARESTPARTICLE);
  buffer.append(pos);
  buffer.broadcast(m_global_rank);

  typedef std::pair<double,int> DistIdPair;
  typedef multimap<int, DistIdPair> RankPairMMap;
  RankPairMMap rankPairMMap;
  m_tml_global_comm.gather(rankPairMMap);

  RankPairMMap::const_iterator it = rankPairMMap.begin();
  DistIdPair nearest = it->second;
  console.Debug()
    << "CLatticeMaster::findParticleNearestTo:"
    << " nearest=" << it->second.second
    << ", distance = " << it->second.first << "\n";

  for (
    it++;
    it != rankPairMMap.end();
    it++
  )
  {
    console.Debug()
      << "CLatticeMaster::findParticleNearestTo:"
      << " particle=" << it->second.second
      << ", distance = " << it->second.first << "\n";

    if (it->second.first < nearest.first)
    {
      console.Debug()
        << "CLatticeMaster::findParticleNearestTo:"
        << " nearest=" << it->second.second
        << ", distance = " << it->second.first << "\n";
      nearest = it->second;
    }
  }
  barrier.wait(
    (std::string("cmd=")+StringUtil::toString(CMD_FINDNEARESTPARTICLE)).c_str()
  );
  console.Debug() << "CLatticeMaster::findParticleNearestTo: exit\n";
  return nearest.second;
}

Vec3 CLatticeMaster::getParticlePosn(int particleId)
{
  console.Debug() << "CLatticeMaster::getParticlePosn: enter\n";
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer buffer(m_global_comm);
  CMPIBarrier barrier(m_global_comm);

  cmd_buffer.broadcast(CMD_GETPARTICLEPOSN);
  buffer.append(particleId);
  buffer.broadcast(m_global_rank);

  typedef std::pair<int,Vec3> IdVec3Pair;
  typedef multimap<int, IdVec3Pair> RankPairMMap;
  RankPairMMap rankPairMMap;
  m_tml_global_comm.gather(rankPairMMap);

  IdVec3Pair idPosnPair(-1,Vec3());
  for (
    RankPairMMap::const_iterator it = rankPairMMap.begin();
    it != rankPairMMap.end();
    it++
  )
  {
    if ((it->second.first) > 0)
    {
      idPosnPair = it->second;
    }
  }
  barrier.wait(
    (std::string("cmd=") + StringUtil::toString(CMD_GETPARTICLEPOSN)).c_str()
  );
  console.Debug() << "CLatticeMaster::getParticlePosn: exit\n";
  return idPosnPair.second;
}

/*!
  Set all particles with a given tag to be non-dynamic, i.e. to have infinite mass

  \param tag the tag
*/
void CLatticeMaster::setParticleNonDynamic(int tag)
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer buffer(m_global_comm);
  CMPIBarrier barrier(m_global_comm);

  cmd_buffer.broadcast(CMD_PSETND);
  buffer.append(tag);
  buffer.broadcast(m_global_rank);
  barrier.wait("setParticleNonDynamic");
}

/*!
  Set all particles with a given tag to be non-rotational, i.e. to have infinite rotational inertia

  \param tag the tag
*/
void CLatticeMaster::setParticleNonRot(int tag)
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer buffer(m_global_comm);
  CMPIBarrier barrier(m_global_comm);

  cmd_buffer.broadcast(CMD_PSETNR);
  buffer.append(tag);
  buffer.broadcast(m_global_rank);
  barrier.wait("setParticleNonRot");
}

/*!
  Set all particles with a given tag to be linear non-dynamic, i.e. switch off velocity updates

  \param tag the tag
*/
void CLatticeMaster::setParticleNonTrans(int tag)
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer buffer(m_global_comm);
  CMPIBarrier barrier(m_global_comm);

  cmd_buffer.broadcast(CMD_PSETNT);
  buffer.append(tag);
  buffer.broadcast(m_global_rank);
  barrier.wait("setParticleNonTrans");
}

/*!
    reset orientation of all particles with a given tag & mask

    \param tag 
    \param mask
*/
void CLatticeMaster::resetParticleRotation(int tag,int mask)
{
    CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
    CVarMPIBuffer buffer(m_global_comm);
    CMPIBarrier barrier(m_global_comm);

    cmd_buffer.broadcast(CMD_RPROT);
    buffer.append(tag);
    buffer.append(mask);
    buffer.broadcast(m_global_rank);
    barrier.wait("resetParticleRotation");
}

/*!
  set the velocity of a particle

  \param id  the id of the particle to be moved
  \param V the velocity
*/
void CLatticeMaster::setParticleVel(int id,const Vec3& V)
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer buffer(m_global_comm);
  CMPIBarrier barrier(m_global_comm);

  cmd_buffer.broadcast(CMD_PVEL);
  buffer.append(id);
  buffer.append(V) ;
  buffer.broadcast(m_global_rank);
  barrier.wait("setParticleVel");
}

/*!
  set the velocity of all particles with a given tag

  \param tag the tag
  \param V the velocity
*/

void CLatticeMaster::setTaggedParticleVel(int tag,const Vec3& V)
{
   CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
   CVarMPIBuffer buffer(m_global_comm);
   CMPIBarrier barrier(m_global_comm);

   cmd_buffer.broadcast(CMD_PTVEL);
   buffer.append(tag);
   buffer.append(V) ;
   buffer.broadcast(m_global_rank);
   barrier.wait("setTaggedParticleVel");
}

/*!
  set the density (i.e. adjust mass) of all particles with a given tag

  \param tag the tag
  \param mask the tag mask
  \param rho the density
*/
void CLatticeMaster::setParticleDensity(int tag,int mask,double rho)
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer buffer(m_global_comm);
  CMPIBarrier barrier(m_global_comm);

  cmd_buffer.broadcast(CMD_PDENS);
  buffer.append(tag);
  buffer.append(mask);
  buffer.append(rho) ;
  buffer.broadcast(m_global_rank);
  barrier.wait("setParticleDensity");
}

/*!
  Call the SubLattice function to set the angular velocity of a particle.
  If the SubLattice is not a RotSubLattice the called function is a NOP.

  \param id the id of the particle
  \param A the angular velocity
*/
void CLatticeMaster::setParticleAngVel(int id,const Vec3& A)
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer buffer(m_global_comm);
  CMPIBarrier barrier(m_global_comm);

  cmd_buffer.broadcast(CMD_PANGVEL);
  buffer.append(id);
  buffer.append(A) ;
  buffer.broadcast(m_global_rank);
  barrier.wait("setParticleAngVel");
}

/*!
  move a wall by given vector

  \param name the name of the wall to be moved
  \param d the movement vector
*/
void CLatticeMaster::moveWallBy(const string& name,const Vec3& d)
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer buffer(m_global_comm);
  CMPIBarrier barrier(m_global_comm);

  cmd_buffer.broadcast(CMD_WMOVE);
  buffer.append(name.c_str());
  buffer.append(d);
  buffer.broadcast(m_global_rank);
  barrier.wait("moveWallBy");
}

/*!
  move a sphere body by given vector

  \param name the name of the sphere body to be moved
  \param d the movement vector
*/
void CLatticeMaster::moveSphereBodyBy(const string& name,const Vec3& d)
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer buffer(m_global_comm);
  CMPIBarrier barrier(m_global_comm);

  cmd_buffer.broadcast(CMD_SPHEREBODYMOVE);
  buffer.append(name.c_str());
  buffer.append(d);
  buffer.broadcast(m_global_rank);
  barrier.wait("moveSphereBodyBy");
}

/*!
  set a wall Normal by given vector

  \param id the nr of the wall normal to change
  \param d the normal
*/
void CLatticeMaster::setWallNormal(const string& name,const Vec3& d)
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer buffer(m_global_comm);
  CMPIBarrier barrier(m_global_comm);

  cmd_buffer.broadcast(CMD_WNORM);
  buffer.append(name.c_str());
  buffer.append(d);
  buffer.broadcast(m_global_rank);
  barrier.wait("setWallNormal");
}

/*!
  apply a given force to a wall

  \param id the nr of the wall
  \param f the force
*/
void CLatticeMaster::applyForceToWall(const string& name,const Vec3& f)
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer buffer(m_global_comm);
  CMPIBarrier barrier(m_global_comm);

  cmd_buffer.broadcast(CMD_WFORCE);
  buffer.append(name.c_str());
  buffer.append(f);
  buffer.broadcast(m_global_rank);
  barrier.wait("applyForceToWall");
}

/*!
  set velocity of a wall. Only affects ViscWall, i.e. wall position
  doesn't get updatet !!

  \param id the nr of the wall
  \param v the velocity
*/
void CLatticeMaster::setVelocityOfWall(const string& name,const Vec3& v)
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer buffer(m_global_comm);
  CMPIBarrier barrier(m_global_comm);

  cmd_buffer.broadcast(CMD_WVEL);
  buffer.append(name.c_str());
  buffer.append(v);
  buffer.broadcast(m_global_rank);
  barrier.wait("setVelocityOfWall");
}

/*!
  perform a single time step
 */
void CLatticeMaster::oneStep()
{
  m_pTimers->start("OneStepForceAndExchange");
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CMPIBarrier barrier(m_global_comm);

  m_pTimers->start("ForceCalculation");
  cmd_buffer.broadcast(CMD_CALC); // local calculations
  barrier.wait("oneStep.1");
  m_pTimers->stop("ForceCalculation");
  
  m_pTimers->start("BoundaryDataExchange");
  cmd_buffer.broadcast(CMD_XCHANGE); // exchange boundary data
  barrier.wait("oneStep.2");
  m_pTimers->stop("BoundaryDataExchange");
  
  //fluid contents
  if(m_fluidinitiated && m_t*m_dt>=m_t_f*m_dt_f){
    m_pTimers->start("UpdateFluid");
    cmd_buffer.broadcast(CMD_UPDATE_FLUID);
    barrier.wait("UpdateFluid");
    m_pTimers->stop("UpdateFluid");  

    m_pTimers->start("ExchangeCells");
    int nt=m_t-m_t_old;
    CVarMPIBuffer buffer(m_global_comm);
    cmd_buffer.broadcast(CMD_EXCHANGE_CELLS);
    buffer.append(nt);
    buffer.broadcast(m_global_rank);
    barrier.wait("ExchangeCells");
    m_pTimers->stop("ExchangeCells");  

    m_pTimers->start("CalculatePressure");
    calculatePressure();
    m_pTimers->stop("CalculatePressure");
    m_t_old=m_t;
    m_t_f++;
  }

  console.Debug() << "=== TIME STEP ===\n";
  //-- times --
  console.Info() << "Force calculation  " << m_pTimers->getTiming("ForceCalculation") << "\n";
  console.Info() << "Boundary exchange  " << m_pTimers->getTiming("BoundaryDataExchange") << "\n";
  console.Info() << "Update fluid  " << m_pTimers->getTiming("UpdateFluid") << "\n"; //fluid contents
  console.Info() << "Calculate Pressure  " << m_pTimers->getTiming("CalculatePressure") << "\n"; //fluid contents
  m_pTimers->stop("OneStepForceAndExchange");
  console.Info() << "one step           " << m_pTimers->getTiming("OneStepForceAndExchange") << "\n";
  //---
}

/****fluid contents: begin****/
bool sortOnZ(const pair<Vec3,double> a, const pair<Vec3,double> b) {
  Vec3 a_index=a.first;
  Vec3 b_index=b.first;
  return a_index.Z() < b_index.Z();
}

bool sortOnY(const pair<Vec3,double> a, const pair<Vec3,double> b) {
  Vec3 a_index=a.first;
  Vec3 b_index=b.first;
  return a_index.Y() < b_index.Y();
}

bool sortOnX(const pair<Vec3,double> a, const pair<Vec3,double> b) {
  Vec3 a_index=a.first;
  Vec3 b_index=b.first;
  return a_index.X() < b_index.X();
}

vector<pair<Vec3,double> > CLatticeMaster::sortCoefficient(vector<pair<Vec3,double> > coefficient) 
{
  vector<pair<Vec3,double> > coeffi = coefficient;
  std::stable_sort(coeffi.begin(), coeffi.end(), sortOnZ);
  std::stable_sort(coeffi.begin(), coeffi.end(), sortOnY);
  std::stable_sort(coeffi.begin(), coeffi.end(), sortOnX);
  return coeffi;
}

void CLatticeMaster::calculatePressure()
{
  //getting coefficients from slaves
  m_pTimers->start("collectcoefficients");
  vector<pair<Vec3,double> > W,E,N,S,D,U,C,B;
  string coeffi_name="c_W";
  W=collectSingleCoeffi(coeffi_name);
  coeffi_name="c_E";
  E=collectSingleCoeffi(coeffi_name);
  coeffi_name="c_N";
  N=collectSingleCoeffi(coeffi_name);
  coeffi_name="c_S";
  S=collectSingleCoeffi(coeffi_name);
  coeffi_name="c_D";
  D=collectSingleCoeffi(coeffi_name);
  coeffi_name="c_U";
  U=collectSingleCoeffi(coeffi_name);
  coeffi_name="c_C";
  C=collectSingleCoeffi(coeffi_name);
  coeffi_name="c_B";
  B=collectSingleCoeffi(coeffi_name);
  m_pTimers->stop("collectcoefficients");
  //sorting coefficients according to cell positions
  m_pTimers->start("sortcoefficients");
  vector<pair<Vec3,double> > sorted_W,sorted_E,sorted_N,sorted_S,sorted_D,sorted_U,sorted_C,sorted_B;
  sorted_W=sortCoefficient(W);
  sorted_E=sortCoefficient(E);
  sorted_N=sortCoefficient(N);
  sorted_S=sortCoefficient(S);
  sorted_D=sortCoefficient(D);
  sorted_U=sortCoefficient(U);
  sorted_C=sortCoefficient(C);
  sorted_B=sortCoefficient(B);

  int minX,minY,minZ,maxX,maxY,maxZ;
  vector<pair<Vec3,double> >::iterator iter_begin=sorted_W.begin();
  Vec3 index_begin=iter_begin->first;
  minX=index_begin.X();minY=index_begin.Y();minZ=index_begin.Z();
  vector<pair<Vec3,double> >::iterator iter_end=sorted_W.end()-1;
  Vec3 index_end=iter_end->first;
  maxX=index_end.X();maxY=index_end.Y();maxZ=index_end.Z();
  int nx=maxX-minX+1;
  int ny=maxY-minY+1;
  int nz=maxZ-minZ+1;
  m_pTimers->stop("sortcoefficients");

  m_pTimers->start("GMRES_solver");
  if(nx*ny*nz%(m_global_size-1)!=0){
    cerr<< "!!Error:Matrix rows("<<nx*ny*nz<<") can not be striped evenly by number of workers("<<m_global_size-1<<") Please choose an appropriate work number!"<<endl;
    abort();
  }
  m_master_solver=new GMRESSolverMaster(&m_global_comm,&m_tml_global_comm,sorted_W,sorted_E,sorted_N,sorted_S,sorted_D,sorted_U,sorted_C,sorted_B,nx,ny,nz);
  vector<pair<Vec3,double> > dis_P;
  dis_P=m_master_solver->MatrixSolving();
  m_pTimers->stop("GMRES_solver");

  m_pTimers->start("distributePressure");
  distributePressure(dis_P);
  m_pTimers->start("distributePressure");

  console.Debug() << "!!!!! FLUID CALCULATION !!!!!\n";
  //-- times --
  console.Info() << "collectcoefficients " << m_pTimers->getTiming("collectcoefficients") << "\n";
  console.Info() << "sortcoefficients " << m_pTimers->getTiming("sortcoefficients") << "\n";
  console.Info() << "GMRES_solver  " << m_pTimers->getTiming("GMRES_solver") << "\n";
  console.Info() << "distributePressure  " << m_pTimers->getTiming("distributePressure") << "\n";
  
  //clean up
  vector<pair<Vec3,double> > zero;
  W.swap(zero);E.swap(zero);N.swap(zero);S.swap(zero);D.swap(zero);U.swap(zero);C.swap(zero);B.swap(zero);
  sorted_W.swap(zero);sorted_E.swap(zero);sorted_N.swap(zero);sorted_S.swap(zero);sorted_D.swap(zero);sorted_U.swap(zero);sorted_C.swap(zero);sorted_B.swap(zero);
  dis_P.clear();
  delete m_master_solver;
}


vector<pair<Vec3,double> > CLatticeMaster::collectSingleCoeffi(const string& coeffi_name)
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CMPIBarrier barrier(m_global_comm);
  cmd_buffer.broadcast(CMD_SEND_COEFFI);

  m_tml_global_comm.broadcast_cont(coeffi_name);
  multimap<int,pair<Vec3,double> > temp_mm;
  m_tml_global_comm.gather(temp_mm);
  barrier.wait((std::string("sendCoefficient_")+coeffi_name).c_str());

  vector<pair<Vec3,double> > coeffi_vector;
  for(multimap<int,pair<Vec3,double> >::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    coeffi_vector.push_back(iter->second);
  };
  return coeffi_vector;
}



void CLatticeMaster::distributePressure(vector<pair<Vec3,double> > P)
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CMPIBarrier barrier(m_global_comm);
  cmd_buffer.broadcast(CMD_RECV_PRESSURE);
  m_tml_global_comm.broadcast_cont(P);
  barrier.wait("Distribute pressure");
}
/****fluid contents: end****/


/*!
  neighbor search. A check if the search is necessary is performed first

  \param force if true, force neighborsearch even if displacement is below threshold
*/
void CLatticeMaster::searchNeighbors(bool force)
{
  console.Debug() << "CLatticeMaster::searchNeighbors\n" ;

  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CMPIBarrier barrier(m_global_comm);

  bool need_search=true;

  if(!force){
    need_search=checkNeighbors();
  }
  if(need_search){
    cmd_buffer.broadcast(CMD_NSEARCH);
    barrier.wait("searchNeighbors");
  }
  console.Debug() << "end CLatticeMaster::searchNeighbors\n" ;
}

/*!
  Test if neighbor search is necessary.

*/
bool CLatticeMaster::checkNeighbors()
{
  console.Debug() << "CLatticeMaster::checkNeighbors()\n";
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CMPIBarrier barrier(m_global_comm);
  CMPISGBufferRoot buffer(m_global_comm,8);

  bool res=false;

  cmd_buffer.broadcast(CMD_CHECKNEIGHBORS);
  buffer.gather();
  int i=1;
  while((i<m_global_size)&&(!res)){
    int b=buffer.pop_int(i);
    res=(b==1); // replace with pop_bool -> to be implemented
    i++;
  }
  barrier.wait("checkNeighbors");
  console.Debug() << "end CLatticeMaster::checkNeighbors(), res= " << res << "\n";

  return res;
}

/*!
 */
void CLatticeMaster::updateInteractions()
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CMPIBarrier barrier(m_global_comm);

  console.XDebug() << "CLatticeMaster::updateInteractions()\n";
  // send the command
  cmd_buffer.broadcast(CMD_UPDATE);

  barrier.wait("updateInteractions");
  console.XDebug() << "end CLatticeMaster::updateInteractions()\n";
}

/*!
  add a scalar particle field to the list of fields to be saved

  \param filename the name of the file the field is saved into
  \param fieldname the name of the field
  \param savetype output file format
  \param t_0 first timestep to be saved
  \param t_end last timestep to be saved
  \param dt timesteps between saves
*/
void CLatticeMaster::addScalarParticleSaveField(
  const string& filename,
  const string& fieldname,
  const string& savetype,
  int t_0,
  int t_end,
  int dt
)
{
  console.Debug()
    << "CLatticeMaster::addScalarParticleSaveField("
    << filename
    << ","
    << fieldname
    << ","
    << t_0
    << ","
    << t_end
    << ","
    << dt
    << ")\n";

  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CMPIBarrier barrier(m_global_comm);

  //send cmd to slave
  cmd_buffer.broadcast(CMD_ADD_SPF);

  AFieldMaster* new_fm =
    new ScalarParticleFieldMaster(
      &m_tml_global_comm,
      fieldname,
      filename,
      savetype,
      t_0,
      t_end,
      dt
    );
  m_save_fields.push_back(new_fm);

  barrier.wait("CLatticeMaster::addScalarParticleSaveField");

  console.Debug() << "end CLatticeMaster::addScalarParticleSaveField()\n";
}

/*!
  add a scalar particle field to the list of fields to be saved

  \param filename the name of the file the field is saved into
  \param fieldname the name of the field
  \param savetype output file format
  \param t_0 first timestep to be saved
  \param t_end last timestep to be saved
  \param dt timesteps between saves
  \param tag tag of the particles for which the field is saved
  \param mask the mask used in tag comparisons
*/
void CLatticeMaster::addTaggedScalarParticleSaveField(const string& filename, const string& fieldname,const string& savetype, int t_0,int t_end,int dt,int tag,int mask)
{
  console.Debug() << "CLatticeMaster::addTaggedScalarParticleSaveField(" << filename << ","<< fieldname << "," << t_0 << "," << t_end << "," << dt << ")\n";

  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CMPIBarrier barrier(m_global_comm);

  //send cmd to slave
  cmd_buffer.broadcast(CMD_ADD_SPF);

  AFieldMaster* new_fm=new ScalarParticleFieldMaster(&m_tml_global_comm,fieldname,filename,savetype,t_0,t_end,dt,tag,mask);
  m_save_fields.push_back(new_fm);

  barrier.wait("CLatticeMaster::addTaggedScalarParticleSaveField");

  console.Debug() << "end CLatticeMaster::addTaggedScalarParticleSaveField()\n";
}

/*!
  save the distribution/histogram of a scalar field on tagged particles

  \param filename the name of the file the field is saved into
  \param fieldname the name of the field
  \param savetype output file format (WINDOW or GLOBAL)
  \param t_0 first timestep to be saved
  \param t_end last timestep to be saved
  \param dt timesteps between adding samples into the distribution
  \param t_snap time between snapshot saves
  \param tag tag of the particles for which the field is saved
  \param mask the mask used in tag comparisons
  \param x_0 minimum data size in distribution
  \param x_max maximum data size in distribution
  \param nx number of bins
*/
void CLatticeMaster::addTaggedScalarParticleDistributionSaver(const string& filename,const string& fieldname,const string& savetype,int t_0,int t_end,int dt,int t_snap,int tag,int mask,double x_0,double x_max,int nx)
{
  console.Debug() << "CLatticeMaster::addTaggedScalarParticleDistributionSaver(" << filename << ","<< fieldname << "," << t_0 << "," << t_end << "," << dt << ")\n";

  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CMPIBarrier barrier(m_global_comm);

  //send cmd to slave
  cmd_buffer.broadcast(CMD_ADD_SPF);

  AFieldMaster* new_fm=new ScalarParticleDistributionMaster(&m_tml_global_comm,fieldname,filename,savetype,t_0,t_end,dt,t_snap,x_0,x_max,nx,tag,mask);
  m_save_fields.push_back(new_fm);

  barrier.wait("CLatticeMaster::addTaggedScalarParticleDistributionSaver");

  console.Debug() << "end CLatticeMaster::addTaggedScalarParticleSaveField()\n";
}

/*!
  add a vector particle field to the list of fields to be saved

  \param filename the name of the file the field is saved into
  \param fieldname the name of the field
  \param savetype output file format
  \param t_0 first timestep to be saved
  \param t_end last timestep to be saved
  \param dt timesteps between saves
*/
void CLatticeMaster::addVectorParticleSaveField(const string& filename,const string& fieldname,const string& savetype,int t_0,int t_end,int dt)
{
  console.Debug() << "CLatticeMaster::addVectorParticleSaveField(" << filename << ","<< fieldname << "," << t_0 << "," << t_end << "," << dt << ")\n";

  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CMPIBarrier barrier(m_global_comm);

  //send cmd to slave
  cmd_buffer.broadcast(CMD_ADD_VPF);

  AFieldMaster* new_fm=new VectorParticleFieldMaster(&m_tml_global_comm,fieldname,filename,savetype,t_0,t_end,dt);
  m_save_fields.push_back(new_fm);

  barrier.wait("CLatticeMaster::addVectorParticleSaveField");

  console.Debug() << "end CLatticeMaster::addVectorParticleSaveField()\n";
}

/*!
  add a vector particle field to the list of fields to be saved

  \param filename the name of the file the field is saved into
  \param fieldname the name of the field
  \param savetype output file format
  \param t_0 first timestep to be saved
  \param t_end last timestep to be saved
  \param dt timesteps between saves
  \param tag tag of the particles for which the field is saved
  \param mask the mask used in tag comparisons
*/
void CLatticeMaster::addTaggedVectorParticleSaveField(const string& filename,const string& fieldname,const string& savetype,int t_0,int t_end,int dt,int tag, int mask)
{
  console.Debug() << "CLatticeMaster::addTaggedVectorParticleSaveField(" << filename << ","<< fieldname << "," << t_0 << "," << t_end << "," << dt << "," << tag << "," << mask << ")\n";

  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CMPIBarrier barrier(m_global_comm);

  //send cmd to slave
  cmd_buffer.broadcast(CMD_ADD_VPF);

  AFieldMaster* new_fm=new VectorParticleFieldMaster(&m_tml_global_comm,fieldname,filename,savetype,t_0,t_end,dt,tag,mask);
  m_save_fields.push_back(new_fm);

  barrier.wait("CLatticeMaster::addTaggedVectorParticleSaveField");

  console.Debug() << "end CLatticeMaster::addTaggedVectorParticleSaveField()\n";
}
/*!
  add a vector particle field to the list of fields to be saved

  \param filename the name of the file the field is saved into
  \param fieldname the name of the field
  \param savetype output file format
  \param t_0 first timestep to be saved
  \param t_end last timestep to be saved
  \param dt timesteps between saves
  \param tprms trigger parameters
*/
void CLatticeMaster::addVectorParticleSaveFieldWT(const string& filename,const string& fieldname,const string& savetype,int t_0,int t_end,int dt,const MaxTrigParams &tprms)
{
  console.Debug() << "CLatticeMaster::addVectorParticleSaveField(" << filename << ","<< fieldname << "," << t_0 << "," << t_end << "," << dt << ")\n";

  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CMPIBarrier barrier(m_global_comm);

  //send cmd to slave
  cmd_buffer.broadcast(CMD_ADD_VPF);

  AFieldMaster* new_fm=new TriggeredVectorParticleFieldMaster(&m_tml_global_comm,fieldname,filename,savetype,t_0,t_end,dt,tprms);
  m_save_fields.push_back(new_fm);

  barrier.wait("CLatticeMaster::addVectorParticleSaveField");

  console.Debug() << "end CLatticeMaster::addVectorParticleSaveField()\n";
}

/*!
  add a vector particle field to the list of fields to be saved

  \param filename the name of the file the field is saved into
  \param fieldname the name of the field
  \param savetype output file format
  \param t_0 first timestep to be saved
  \param t_end last timestep to be saved
  \param dt timesteps between saves
  \param tag tag of the particles for which the field is saved
  \param mask the mask used in tag comparisons
  \param tprms trigger parameters
*/
void CLatticeMaster::addTaggedVectorParticleSaveFieldWT(const string& filename,const string& fieldname,const string& savetype,int t_0,int t_end,int dt,int tag, int mask,const MaxTrigParams &tprms)
{
  console.Debug() << "CLatticeMaster::addTaggedVectorParticleSaveField(" << filename << ","<< fieldname << "," << t_0 << "," << t_end << "," << dt << "," << tag << "," << mask << ")\n";

  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CMPIBarrier barrier(m_global_comm);

  //send cmd to slave
  cmd_buffer.broadcast(CMD_ADD_VPF);

  AFieldMaster* new_fm=new TriggeredVectorParticleFieldMaster(&m_tml_global_comm,fieldname,filename,savetype,t_0,t_end,dt,tag,mask,tprms);
  m_save_fields.push_back(new_fm);

  barrier.wait("CLatticeMaster::addTaggedVectorParticleSaveField");

  console.Debug() << "end CLatticeMaster::addTaggedVectorParticleSaveField()\n";
}

/*!
  add a vector wall field to the list of fields to be saved

  \param filename the name of the file the field is saved into
  \param fieldname the name of the field
  \param walls names of the walls
  \param t_0 first timestep to be saved
  \param t_end last timestep to be saved
  \param dt timesteps between saves
*/
void CLatticeMaster::addVectorWallField(const string& filename,
                    const string& fieldname,
                    vector<string> walls,
                    const string& savetype,
                    int t_0,
                    int t_end,
                    int dt)
{
  console.Debug() << "CLatticeMaster::addVectorWallField(" << filename << ","<< fieldname << "," << t_0 << "," << t_end << "," << dt << ")\n";

  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);

  //send command to slaves
  cmd_buffer.broadcast(CMD_ADD_VWF);

  // create master
  VectorWallFieldMaster* new_fm=new VectorWallFieldMaster(&m_tml_global_comm,fieldname,filename,walls,savetype,t_0,t_end,dt);
  m_save_fields.push_back(new_fm);

  m_tml_global_comm.barrier();

  console.Debug() << "end CLatticeMaster::addVectorWallField()\n";
}

/*!
  Setup parameters for restart checkpoints

 \param fileNamePrefix Path prefix for checkpoint files. Multiple snapshot
  files may be generated for a single timestep snapshot.
  \param beginTime      Time to begin checkpointing. Time of first snapshot.
  \param endTime        End time for checkpointing. Time of last snapshot.
  \param timeInterval   Time interval between snapshot file generation.
  \param saveBinary     Saves the data in binary format if true, ascii if false
*/
void CLatticeMaster::performCheckPoints(
  const string &fileNamePrefix,
  int beginTime,
  int endTime,
  int timeInterval,
  int precision
)
{
  m_pCheckPointController->setGeometryInfo(m_geo_info);
  m_pCheckPointController->setCheckPointParams(
    fileNamePrefix,
    beginTime,
    endTime,
    timeInterval,
    false,
    precision
  );
  m_pCheckPointController->setMpiComm(m_global_comm);
}

/*!
  Setup parameters for restart checkpoints written though the master process

 \param fileNamePrefix Path prefix for checkpoint files. Multiple snapshot
  files may be generated for a single timestep snapshot.
  \param beginTime      Time to begin checkpointing. Time of first snapshot.
  \param endTime        End time for checkpointing. Time of last snapshot.
  \param timeInterval   Time interval between snapshot file generation.
  \param saveBinary     Saves the data in binary format if true, ascii if false
*/
void CLatticeMaster::performCheckPointsThroughMaster(
  const string &fileNamePrefix,
  int beginTime,
  int endTime,
  int timeInterval,
  int precision
)
{
  // check if spatial extent of the model has been set and warn if not
  if(!m_bbx_has_been_set){
    console.Warning() << "Setting up Checkpointer : Model Spatial Domain has not been set - checkpoint headers will contain invalid geometry info \n";
  }
  m_pSnapShotController->setGeometryInfo(m_geo_info);
  m_pCheckPointController->setCheckPointParams(
    fileNamePrefix,
    beginTime,
    endTime,
    timeInterval,
    true,
    precision
  );
  m_pCheckPointController->setMpiComm(m_global_comm);
}

/*!
  Initialises parameters for performing model snapshots.

  \param fileNamePrefix Path prefix for checkpoint files. Multiple snapshot
  files may be generated for a single timestep snapshot.
  \param beginTime      Time to begin checkpointing. Time of first snapshot.
  \param endTime        End time for checkpointing. Time of last snapshot.
  \param timeInterval   Time interval between snapshot file generation.
*/
void CLatticeMaster::initSnapShotController(
  const string &fileNamePrefix,
  int beginTime,
  int endTime,
  int timeInterval
)
{
  console.Debug() << "CLatticeMaster::initSnapShotController\n";
  // check if spatial extent of the model has been set and warn if not
  if(!m_bbx_has_been_set){
    console.Warning() << "Setting up Snapshot Controller : Model Spatial Domain has not been set - snapshot headers will contain invalid geometry info \n";
  }
  m_pSnapShotController->setGeometryInfo(m_geo_info);
  m_pSnapShotController->setCheckPointParams(
    fileNamePrefix,
    beginTime,
    endTime,
    timeInterval,
    false
  );
  m_pSnapShotController->setMpiComm(m_global_comm);
  console.Debug() << "end CLatticeMaster::initSnapShotController\n";
}

/*!
  add a scalar interaction field to the list of fields to be saved

  \param filename the name of the file the field is saved into
  \param fieldname the name of the field
  \param igname the name of the interaction group from which the field is taken
  \param savetype output file format
  \param t_0 first timestep to be saved
  \param t_end last timestep to be saved
  \param dt timesteps between saves
  \param checked choose between normal and checked field (defaults to false)
*/
void CLatticeMaster::addScalarInteractionSaveField(
  const string& filename,
  const string& fieldname,
  const string& igtype,
  const string& igname,
  const string& savetype,
  int t_0,
  int t_end,
  int dt,
  bool checked
)
{
  console.Debug()
     << "CLatticeMaster::addScalarInteractionSaveField("
     << filename
     << ","
     << fieldname
     << ","
     << igname
     << ","
     << t_0
     << ","
     << t_end
     << ","
     << dt
     << ")\n";

  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CMPIBarrier barrier(m_global_comm);

  //send cmd to slave
  cmd_buffer.broadcast(CMD_ADD_SIF);

  AFieldMaster* new_fm =
    new ScalarInteractionFieldMaster(
      &m_tml_global_comm,
      fieldname,
      igtype,
      igname,
      filename,
      savetype,
      t_0,
      t_end,
      dt,
      checked
    );
  m_save_fields.push_back(new_fm);

  barrier.wait("CLatticeMaster::addScalarInteractionSaveField");

  console.Debug() << "end CLatticeMaster::addScalarInteractionSaveField()\n";
}

/*!
  add a scalar interaction field to the list of fields to be saved

  \param filename the name of the file the field is saved into
  \param fieldname the name of the field
  \param igname the name of the interaction group from which the field is taken
  \param savetype output file format
  \param t_0 first timestep to be saved
  \param t_end last timestep to be saved
  \param dt timesteps between saves
  \param checked choose between normal and checked field (defaults to false)
*/
void CLatticeMaster::addScalarHistoryInteractionSaveField(const string& filename,const string& fieldname,
    const string& igtype,const string& igname,const string& savetype,int t_0,int t_end,int dt)
{
  console.Debug() << "CLatticeMaster::addScalarHistoryInteractionSaveField(" << filename << "," << fieldname << "," << igname << ","
    << t_0 << "," << t_end << "," << dt << ")\n";

  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CMPIBarrier barrier(m_global_comm);

  //send cmd to slave
  cmd_buffer.broadcast(CMD_ADD_HIF);

  //AFieldMaster* new_fm = new ScalarHistoryInteractionFieldMaster(&m_tml_global_comm,fieldname,igtype,igname,
  //    filename,savetype,t_0,t_end,dt,checked);
  //m_save_fields.push_back(new_fm);

  barrier.wait("CLatticeMaster::addScalarHistoryInteractionSaveField");

  console.Debug() << "end CLatticeMaster::addScalarHistoryInteractionSaveField()\n";
}

/*!
  add a vector field on the triangles of a given trimesh to the list of fields to be saved

  \param filename the name of the file the field is saved into
  \param fieldname the name of the field
  \param meshname the name of the mesh from which the field is taken
  \param savetype the format in which the data is to be saved
  \param t_0 first timestep to be saved
  \param t_end last timestep to be saved
  \param dt timesteps between saves
*/
void CLatticeMaster::addVectorTriangleSaveField(const string& filename,
                        const string& fieldname,
                        const string& meshname,
                        const string& savetype,
                        int t_0,int t_end,int dt)
{
  console.Debug() << "CLatticeMaster::addVectorTriangleSaveField\n";

  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CMPIBarrier barrier(m_global_comm);

  //send cmd to slave
  cmd_buffer.broadcast(CMD_ADD_VTF);

  AFieldMaster* new_fm = new VectorTriangleFieldMaster(&m_tml_global_comm,fieldname,meshname,filename,savetype,t_0,t_end,dt);
  m_save_fields.push_back(new_fm);

  barrier.wait("CLatticeMaster::addVectorTriangleSaveField");

  console.Debug() << "end CLatticeMaster::addVectorTriangleSaveField\n";
}

/*!
  add a scalar field on the triangles of a given trimesh to the list of fields to be saved

  \param filename the name of the file the field is saved into
  \param fieldname the name of the field
  \param meshname the name of the mesh from which the field is taken
  \param savetype the format in which the data is to be saved
  \param t_0 first timestep to be saved
  \param t_end last timestep to be saved
  \param dt timesteps between saves
*/
void CLatticeMaster::addScalarTriangleSaveField(const string& filename,
                        const string& fieldname,
                        const string& meshname,
                        const string& savetype,
                        int t_0,int t_end,int dt)
{
  console.Debug() << "CLatticeMaster::addScalarTriangleSaveField\n";

  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CMPIBarrier barrier(m_global_comm);

  //send cmd to slave
  cmd_buffer.broadcast(CMD_ADD_STF);

  AFieldMaster* new_fm = new ScalarTriangleFieldMaster(&m_tml_global_comm,fieldname,meshname,filename,savetype,t_0,t_end,dt);
  m_save_fields.push_back(new_fm);

  barrier.wait("CLatticeMaster::addScalarTriangleSaveField");

  console.Debug() << "end CLatticeMaster::addScalarTriangleSaveField\n";
}

/*!
  add a vector interaction field to the list of fields to be saved

  \param filename the name of the file the field is saved into
  \param fieldname the name of the field
  \param igname the name of the interaction group from which the field is taken
  \param savetype output file format
  \param t_0 first timestep to be saved
  \param t_end last timestep to be saved
  \param dt timesteps between saves
  \param checked choose between normal and checked field (defaults to false)
*/
void CLatticeMaster::addVectorInteractionSaveField(const string& filename,const string& fieldname,const string& igtype,
                           const string& igname,const string& savetype,int t_0,int t_end,int dt,bool checked)
{
  console.Debug() << "CLatticeMaster::addVectorInteractionSaveField()\n";

  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CMPIBarrier barrier(m_global_comm);

  //send cmd to slave
  cmd_buffer.broadcast(CMD_ADD_VIF);

  AFieldMaster* new_fm = new VectorInteractionFieldMaster(&m_tml_global_comm,fieldname,igtype,igname,filename,savetype,t_0,t_end,dt,checked);
  m_save_fields.push_back(new_fm);

  barrier.wait("CLatticeMaster::addVectorInteractionSaveField");

  console.Debug() << "end CLatticeMaster::addVectorInteractionSaveField()\n";
}


/*!
  add a scalar interaction field to the list of fields to be saved

  \param filename the name of the file the field is saved into
  \param fieldname the name of the field
  \param igname the name of the interaction group from which the field is taken
  \param savetype output file format
  \param t_0 first timestep to be saved
  \param t_end last timestep to be saved
  \param dt timesteps between saves
  \param tag the particle tag
  \param mask the mask used for tag comparisons
    \param checked choice between "full" and "checked" fields
*/
void CLatticeMaster::addTaggedScalarInteractionSaveField(const string& filename, const string& fieldname,const string& igtype,const string& igname,const string& savetype, int t_0,int t_end,int dt,int tag,int mask, bool checked)
{
  console.Debug() << "CLatticeMaster::addTaggedScalarInteractionSaveField(" << filename << ","<< fieldname << "," << igname << "," << t_0 << "," << t_end << "," << dt << ")\n";

  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CMPIBarrier barrier(m_global_comm);

  //send cmd to slave
  cmd_buffer.broadcast(CMD_ADD_SIF);

  AFieldMaster* new_fm=new ScalarInteractionFieldMaster(&m_tml_global_comm,fieldname,igtype,igname,filename,savetype,t_0,t_end,dt,tag,mask,checked);
  m_save_fields.push_back(new_fm);

  barrier.wait("CLatticeMaster::addTaggedScalarInteractionSaveField");

  console.Debug() << "end CLatticeMaster::addTaggedScalarInteractionSaveField()\n";
}

/*!
  Initialisation to run the simulation
*/
void CLatticeMaster::runInit()
{
  if (!m_isInitialized)
  {
    m_pTimers->start("TotalRunTime");

    console.Debug() << "timer resolution : " << MPI_Wtick() << "\n" ;
    m_isInitialized = true ;
    //m_t=0 ;
    m_pTimers->clear();
  }
}

/*!
  Finalize after running the simulation
*/
void CLatticeMaster::runEnd()
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CMPIBarrier barrier(m_global_comm);

  console.Debug() << "total time: " << m_pTimers->getTiming("TotalRunTime") << " seconds for " << m_max_ts-1 << " time steps\n" ;
  cmd_buffer.broadcast(CMD_FINISH);
  barrier.wait("runEnd");
  barrier.wait("runEnd2");
}

void CLatticeMaster::saveTimingData()
{
  if (!getTimingFileName().empty()) {

    m_pTimers->start("WorkerTimingDataSave");
    CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
    CMPIBarrier barrier(m_global_comm);
    cmd_buffer.broadcast(CMD_SAVETIMINGDATA);
    barrier.wait("saveTimingData");
    m_pTimers->stop("WorkerTimingDataSave");

    m_pTimers->appendData(getTimingFileName());
  }
}

void CLatticeMaster::addPreTimeStepRunnable(Runnable &runnable)
{
  m_preRunnableVector.push_back(&runnable);
}

void CLatticeMaster::addPostTimeStepRunnable(Runnable &runnable)
{
  m_postRunnableVector.push_back(&runnable);
}

void CLatticeMaster::runRunnables(RunnableVector::iterator begin, RunnableVector::iterator end)
{
  for (RunnableVector::iterator it = begin; it != end; it++)
  {
    (*it)->run();
  }
}

void CLatticeMaster::runPreRunnables()
{
  runRunnables(
    m_preRunnableVector.begin(),
    m_preRunnableVector.end()
  );
}

void CLatticeMaster::runPostRunnables()
{
  runRunnables(
    m_postRunnableVector.begin(),
    m_postRunnableVector.end()
  );
}

/**
 *  Perform a single time step of the simulation.
 */
void CLatticeMaster::runOneStep()
{
  runInit();

  m_pTimers->start("RunOneStep");

  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CMPIBarrier barrier(m_global_comm);

  // Field saving of initial state (if required)
  if(m_first_time){
    m_pTimers->start("FieldMasterCollect");
    m_pTimers->pause("FieldMasterCollect");
    m_pTimers->start("FieldMasterWrite");
    m_pTimers->pause("FieldMasterWrite");
    for(vector<AFieldMaster*>::iterator iter=m_save_fields.begin();
        iter!=m_save_fields.end();
        iter++){
      if((*iter)->needSave(0)){
        m_pTimers->resume("FieldMasterCollect");
        cmd_buffer.broadcast(CMD_SEND_FIELDS);
        (*iter)->collect();
        barrier.wait("runOneStep.collect");
        m_pTimers->pause("FieldMasterCollect");
        m_pTimers->resume("FieldMasterWrite");
        (*iter)->write();
        m_pTimers->pause("FieldMasterWrite");
      }
    }
    m_pTimers->stop("FieldMasterCollect");
    m_pTimers->stop("FieldMasterWrite");
  }
  // Save a check-point of the initial state (if required)
  m_pTimers->start("NeighbourSearch");
  if(m_first_time){
    searchNeighbors(true);
    m_pCheckPointController->performCheckPoint(0);
    m_pSnapShotController->performSnapShot(0);
    m_first_time=false;
  }
  m_pTimers->stop("NeighbourSearch");

  m_t++ ;
  console.Info() << "Begining time step: " << m_t << "\n" ;
  m_pTimers->start("LoadingMechanism");
  runPreRunnables();
  m_pTimers->stop("LoadingMechanism");

  m_pTimers->start("NeighbourSearch");
  searchNeighbors(false);
  m_pTimers->stop("NeighbourSearch");

  //  barrier.wait("runOneStep.2");
  m_pTimers->start("UpdateInteractions");
  updateInteractions();
  m_pTimers->stop("UpdateInteractions");

  m_pTimers->start("CalcForcesAndExchange");
  oneStep();
  m_pTimers->stop("CalcForcesAndExchange");
  m_pTimers->start("SavingData");
  // saving fields
  /*
  bool app;
  if(m_t==1){
    app=false; // if 1st time step -> new file
  } else {
    app=true; // else append
  }
  */

  // new field saving
  m_pTimers->start("FieldMasterCollect");
  m_pTimers->pause("FieldMasterCollect");
  m_pTimers->start("FieldMasterWrite");
  m_pTimers->pause("FieldMasterWrite");
  for(vector<AFieldMaster*>::iterator iter=m_save_fields.begin();
      iter!=m_save_fields.end();
      iter++){
    if((*iter)->needSave(m_t)){
      m_pTimers->resume("FieldMasterCollect");
      cmd_buffer.broadcast(CMD_SEND_FIELDS);
      (*iter)->collect();
      barrier.wait("runOneStep.collect");
      m_pTimers->pause("FieldMasterCollect");
      m_pTimers->resume("FieldMasterWrite");
      (*iter)->write();
      m_pTimers->pause("FieldMasterWrite");
    }
  }
  m_pTimers->stop("FieldMasterCollect");
  m_pTimers->stop("FieldMasterWrite");

  // Save a check-point (if required)
  m_pTimers->start("CheckPoint");
  // force ntable rebuild before checkpoint
  if(m_pCheckPointController->isCheckPoint(m_t)){
    searchNeighbors(true);
  }
  m_pCheckPointController->performCheckPoint(m_t);
  m_pTimers->stop("CheckPoint");

  // write a snapshot (if required)
  m_pTimers->start("SnapShot");
  m_pSnapShotController->performSnapShot(m_t);
  m_pTimers->stop("SnapShot");

  m_pTimers->stop("SavingData");

  runPostRunnables();

  m_pTimers->stop("RunOneStep");

  m_pTimers->stop("TotalRunTime", true);

  saveTimingData();

  console.Info() << "updateInteractions : " << m_pTimers->getTiming("UpdateInteractions") << " seconds\n" ;
  console.Info() << "searchNeighbors    : " <<  m_pTimers->getTiming("NeighbourSearch") << " seconds\n" ;
  console.Info() << "one step           : " <<  m_pTimers->getTiming("CalcForcesAndExchange") << " seconds\n" ;
  console.Info() << "saving data        : " <<  m_pTimers->getTiming("SavingData") << " seconds\n" ;

  console.Info() << "End of time step: " << m_t << "\n" ;
}

/*!
  run the simulation
*/
void CLatticeMaster::run()
{
  while (m_t < m_max_ts) {

    if((int(m_max_ts/10)!=0 && m_t%(int(m_max_ts/10))==0))
      cout<<"----TIME STEP: "<<m_t<<"/"<<m_max_ts<<"----"<<endl;

    runOneStep();
  }
  runEnd();
}
/*
CLatticeMaster::ParticleIdPairVector
CLatticeMaster::getBondGroupIdPairs(
  const std::string &groupName
)
{
  console.XDebug() << "CLatticeMaster::getBondGroupIdPairs(): enter\n";
  BroadcastCommand cmd(getGlobalRankAndComm(), CMD_GETBONDGROUPIDPAIRS);
  cmd.append(groupName.c_str());
  cmd.broadcastCommand();
  cmd.broadcastBuffer();
  typedef std::pair<double,int> DistIdPair;
  typedef multimap<int, ParticleIdPair> ParticleIdPairMMap;
  ParticleIdPairMMap idPairMMap;
  m_tml_global_comm.gather(idPairMMap);
  ParticleIdPairVector idPairVector;
  idPairVector.reserve(idPairMMap.size());
  std::transform(
    idPairMMap.begin(),
    idPairMMap.end(),
    std::back_insert_iterator<ParticleIdPairVector>(idPairVector),
    ext::select2nd<ParticleIdPairMMap::value_type>()
  );
  cmd.wait("CMD_GETBONDGROUPIDPAIRS");
  console.XDebug() << "CLatticeMaster::getBondGroupIdPairs(): exit\n";
  return idPairVector;
}
*/

/*!
  Create and add a new bonded IG

  \param name the name of the bonded IG
  \param k spring constant of the created bonds
  \param dist maximum distance between two particles to create bond
  \param break_dist (relative) breaking distance of bonds

*/
void CLatticeMaster::addBondedIG(const CBondedIGP &prms)
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer pbuffer(m_global_comm,20);
  CMPIBarrier barrier(m_global_comm);

  console.XDebug() << "CLatticeMaster::addBondedIG()\n ";

  // send the command
  cmd_buffer.broadcast(CMD_ADDBONDEDIG);
  // send parameters
  pbuffer.append(prms.tag);
  pbuffer.append(prms.getName().c_str());
  pbuffer.append(prms.k);
  pbuffer.append(prms.rbreak);
  pbuffer.append(static_cast<int>(prms.m_scaling));

  pbuffer.broadcast(0);

  barrier.wait("addBondedIG()");
  console.XDebug() << "end CLatticeMaster::addBondedIG()" << "\n";
}

/*!
  Create and add a new bonded IG with force limit

  \param name the name of the bonded IG
  \param k spring constant of the created bonds
  \param dist maximum distance between two particles to create bond
  \param break_dist (relative) breaking distance of bonds
  \param maxforce the maximum force
*/
void CLatticeMaster::addCappedBondedIG(int tag,const string& name,double k,double break_dist,double maxforce)
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer pbuffer(m_global_comm,20);
  CMPIBarrier barrier(m_global_comm);

  console.XDebug() << "CLatticeMaster::addCappedBondedIG()\n ";

  // send the command
  cmd_buffer.broadcast(CMD_ADDCAPPEDBONDEDIG);
  // send parameters
  pbuffer.append(tag);
  pbuffer.append(name.c_str());
  pbuffer.append(k);
  pbuffer.append(break_dist);
  pbuffer.append(maxforce);
  pbuffer.broadcast(0);

  // broadcast connection data
  // m_tml_global_comm.broadcast_cont(m_temp_conn[tag]);

  barrier.wait("addCappedBondedIG()");
  console.XDebug() << "end CLatticeMaster::addCappedBondedIG()" << "\n";
}

/*!
  Create and add a new short bonded IG

  \param name the name of the bonded IG
  \param k spring constant of the created bonds
  \param dist maximum distance between two particles to create bond
  \param break_dist (relative) breaking distance of bonds

*/
void CLatticeMaster::addShortBondedIG(int tag,const string& name,double k,double break_dist)
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer pbuffer(m_global_comm,20);
  CMPIBarrier barrier(m_global_comm);

  console.XDebug() << "CLatticeMaster::addShortBondedIG()\n ";

  // send the command
  cmd_buffer.broadcast(CMD_ADDSHORTBONDEDIG);
  // send parameters
  pbuffer.append(tag);
  pbuffer.append(name.c_str());
  pbuffer.append(k);
  pbuffer.append(break_dist);

  pbuffer.broadcast(0);

  // broadcast connection data

  //console.XDebug() << "broadcast " << m_temp_conn[0].size() << " connection data" << "\n";
  //m_tml_global_comm.broadcast_cont(m_temp_conn[tag]);
  barrier.wait("addShortBondedIG()");
  console.XDebug() << "end CLatticeMaster::addShortBondedIG()" << "\n";
}

void CLatticeMaster::addRotBondedIG(
  int    tag,
  const  string& name,
  double kr,
  double ks,
  double kt,
  double kb,
  double max_nForce,
  double max_shForce,
  double max_tMoment,
  double max_bMoment,
  bool   scaling,
  bool   meanR_scaling,
  double truncated
)
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer pbuffer(m_global_comm, 256);  // 20 need change ?
  CMPIBarrier barrier(m_global_comm);

  console.XDebug() << "CLatticeMaster::addRotBondedIG()\n ";

  // send the command
  cmd_buffer.broadcast(CMD_ADDROTBONDEDIG);
  // send parameters
  pbuffer.append(tag);
  pbuffer.append(name.c_str());
  pbuffer.append(kr);
  pbuffer.append(ks);
  pbuffer.append(kt);
  pbuffer.append(kb);
  pbuffer.append(max_nForce);
  pbuffer.append(max_shForce);
  pbuffer.append(max_tMoment);
  pbuffer.append(max_bMoment);
  pbuffer.append(static_cast<int>(scaling));
  pbuffer.append(static_cast<int>(meanR_scaling));
  pbuffer.append(truncated);

  pbuffer.broadcast(0);

  // broadcast connection data
  //m_tml_global_comm.broadcast_cont(m_temp_conn[tag]);

  barrier.wait("addRotBondedIG()");
  console.XDebug() << "end CLatticeMaster::addRotBondedIG()" << "\n";
}

/*!
    Add brittle beam interaction group with stress failure criterion
*/
void CLatticeMaster::addBrittleBeamSCIG(const BrittleBeamSCIGP& prms)
{ 
    CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
    CVarMPIBuffer pbuffer(m_global_comm, 256);  // 20 need change ?
    CMPIBarrier barrier(m_global_comm);

    //console.XDebug() << "CLatticeMaster::addBrittleBeamSCIG()\n ";
    std::cout << "CLatticeMaster::addBrittleBeamSCIG()\n ";
    
    // send the command
    cmd_buffer.broadcast(CMD_ADDBRITTLEBEAMSCIG);
    // send parameters
    pbuffer.append(prms.tag);
    pbuffer.append(prms.getName().c_str());
    pbuffer.append(prms.kr);
    pbuffer.append(prms.ks);
    pbuffer.append(prms.kt);
    pbuffer.append(prms.kb);
    pbuffer.append(prms.cohesion);
    pbuffer.append(prms.tCutoff);
    pbuffer.append(prms.fAngle);
    pbuffer.append(static_cast<int>(prms.scaling));

    pbuffer.broadcast(0);

    barrier.wait("addRotBondedIG()");
    //console.XDebug() << "end CLatticeMaster::addBrittleBeamSCIG()" << "\n";
    std::cout << "end CLatticeMaster::addBrittleBeamSCIG()" << "\n";
}

/*!
    Add brittle beam interaction group with Ding-Zhang failure criterion
*/
void CLatticeMaster::addBrittleBeamDZCIG(const BrittleBeamDZCIGP& prms)
{
    CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
    CVarMPIBuffer pbuffer(m_global_comm, 256);  // 20 need change ?
    CMPIBarrier barrier(m_global_comm);

    //console.XDebug() << "CLatticeMaster::addBrittleBeamDZCIG()\n ";
    std::cout << "CLatticeMaster::addBrittleBeamDZCIG()\n ";

    // send the command
    cmd_buffer.broadcast(CMD_ADDBRITTLEBEAMDZCIG);
    // send parameters
    pbuffer.append(prms.tag);
    pbuffer.append(prms.getName().c_str());
    pbuffer.append(prms.kr);
    pbuffer.append(prms.ks);
    pbuffer.append(prms.kt);
    pbuffer.append(prms.kb);
    pbuffer.append(prms.cohesion);
    pbuffer.append(prms.tCutoff);
    pbuffer.append(prms.cCutoff);
    pbuffer.append(prms.fAngle);
    pbuffer.append(prms.beta1);
    pbuffer.append(prms.beta2);
    pbuffer.append(static_cast<int>(prms.scaling));

    pbuffer.broadcast(0);

    barrier.wait("addRotBondedIG()");
    //console.XDebug() << "end CLatticeMaster::addBrittleBeamDZCIG()" << "\n";
    std::cout << "end CLatticeMaster::addBrittleBeamDZCIG()" << "\n";
}

void CLatticeMaster::addRotThermBondedIG(
  const CRotThermBondedIGP &prms
)
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer pbuffer(m_global_comm, 256);  // 20 need change ?
  CMPIBarrier barrier(m_global_comm);

  console.XDebug() << "CLatticeMaster::addRotThermBondedIG()\n ";

  // send the command
  cmd_buffer.broadcast(CMD_ADDROTTHERMBONDEDIG);
  // send parameters
  pbuffer.append(prms.tag);
  pbuffer.append(prms.getName().c_str());
  pbuffer.append(prms.kr);
  pbuffer.append(prms.ks);
  pbuffer.append(prms.kt);
  pbuffer.append(prms.kb);
  pbuffer.append(prms.max_nForce);
  pbuffer.append(prms.max_shForce);
  pbuffer.append(prms.max_tMoment);
  pbuffer.append(prms.max_bMoment);
  pbuffer.append(prms.diffusivity);

  pbuffer.broadcast(0);

  // broadcast connection data
  //m_tml_global_comm.broadcast_cont(m_temp_conn[tag]);

  barrier.wait("addRotThermBondedIG()");
  console.XDebug() << "end CLatticeMaster::addRotThermBondedIG()" << "\n";
}

int CLatticeMaster::getNumParticles()
{
  BroadcastCommand cmd(getGlobalRankAndComm(), CMD_GETNUMPARTICLES);
  cmd.broadcastCommand();

  // get slave dimensions and coordinates
  typedef multimap<int,int> RankNumParticlesMMap;
  RankNumParticlesMMap numParticlesMMap;
  m_tml_global_comm.gather(numParticlesMMap);

  int numParticles = 0;
  for (
    RankNumParticlesMMap::const_iterator it = numParticlesMMap.begin();
    it != numParticlesMMap.end();
    it++
  )
  {
    numParticles += it->second;
  }
  cmd.wait("CLatticeMaster::getNumParticles");
  return numParticles;
}

class IGPCommand : public BroadcastCommand
{
public:
  IGPCommand(const MpiRankAndComm &globalRankAndComm)
    : BroadcastCommand(globalRankAndComm, CMD_ADDPIG)
  {
  }

  IGPCommand(const MpiRankAndComm &globalRankAndComm, int commandId)
    : BroadcastCommand(globalRankAndComm, commandId)
  {
  }

  void appendIGP(const AIGParam &prms)
  {
    appendTypeAndName(prms);
  }

private:
};

void CLatticeMaster::addPairIG(const CElasticIGP &prms)
{
  IGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.m_k);
  igpCmd.append(static_cast<int>(prms.m_scaling));
  igpCmd.broadcast();
}

void CLatticeMaster::addPairIG(const CFrictionIGP &prms)
{
  IGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.k);
  igpCmd.append(prms.mu);
  igpCmd.append(prms.k_s);
  igpCmd.append(prms.dt);
  igpCmd.append(static_cast<int>(prms.m_scaling));
  igpCmd.broadcast();
}

void CLatticeMaster::addPairIG(const CSpringDashpotFrictionIGP &prms)
{
  IGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.youngModulus);
  igpCmd.append(prms.poissonRatio);
  igpCmd.append(prms.cor);
  igpCmd.append(prms.mu);
  igpCmd.append(prms.dt);
  igpCmd.broadcast();
}

void CLatticeMaster::addPairIG(const FractalFrictionIGP &prms)
{
  IGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.k);
  igpCmd.append(prms.mu_0);
  igpCmd.append(prms.k_s);
  igpCmd.append(prms.dt);
  igpCmd.append(prms.x0);
  igpCmd.append(prms.y0);
  igpCmd.append(prms.dx);
  igpCmd.append(prms.dy);
  igpCmd.append(prms.nx);
  igpCmd.append(prms.ny);
  for (int i = 0; i < (prms.nx*prms.ny); i++) {
    igpCmd.append((prms.mu.get())[i]);
  }
  igpCmd.broadcast();
}

void CLatticeMaster::addPairIG(const CAdhesiveFrictionIGP &prms)
{
  IGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.k);
  igpCmd.append(prms.mu);
  igpCmd.append(prms.k_s);
  igpCmd.append(prms.dt);
  igpCmd.append(prms.r_cut);
  igpCmd.broadcast();
}

void CLatticeMaster::addPairIG(const CLinearDashpotIGP & prms)
{
  IGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.m_damp);
  igpCmd.append(prms.m_cutoff);
  igpCmd.broadcast();
}

void CLatticeMaster::addPairIG(const CHertzianElasticIGP &prms)
{
  IGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.m_E);
  igpCmd.append(prms.m_nu);
  igpCmd.broadcast();
}

void CLatticeMaster::addPairIG(const CHertzianViscoElasticFrictionIGP &prms)
{
  IGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.m_A);
  igpCmd.append(prms.m_E);
  igpCmd.append(prms.m_nu);
  igpCmd.append(prms.mu);
  igpCmd.append(prms.k_s);
  igpCmd.append(prms.dt);
  igpCmd.broadcast();
}

void CLatticeMaster::addPairIG(const CHertzianViscoElasticIGP &prms)
{
  IGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.m_A);
  igpCmd.append(prms.m_E);
  igpCmd.append(prms.m_nu);
  igpCmd.broadcast();
}

void CLatticeMaster::addPairIG(const CHertzMindlinIGP &prms)
{
  IGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.m_E);
  igpCmd.append(prms.m_nu);
  igpCmd.append(prms.mu);
  igpCmd.append(prms.dt);
  igpCmd.broadcast();
}

void CLatticeMaster::addPairIG(const CHertzMindlinViscoIGP &prms)
{
  IGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.m_E);
  igpCmd.append(prms.m_nu);
  igpCmd.append(prms.mu);
  igpCmd.append(prms.m_COR);
  igpCmd.append(prms.dt);
  igpCmd.broadcast();
}

//---
// rotational interactions
//---

class RotIGPCommand : public IGPCommand
{
public:
  RotIGPCommand(const MpiRankAndComm &globalRankAndComm) : IGPCommand(globalRankAndComm, CMD_ADDPIG)
  {
  }

  void appendIGP(const AIGParam &prms)
  {
    appendTypeAndName(prms);
  }

private:
};

void CLatticeMaster::addPairIG(const CRotElasticIGP &prms)
{
  RotIGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.m_kr);
  igpCmd.append(static_cast<int>(prms.m_scaling));
  igpCmd.broadcast();
}

void CLatticeMaster::addPairIG(const CRotFrictionIGP &prms)
{
  RotIGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.k);
  igpCmd.append(prms.mu_s);
  igpCmd.append(prms.mu_d);
  igpCmd.append(prms.k_s);
  igpCmd.append(prms.dt);
  igpCmd.append(static_cast<int>(prms.scaling));
  igpCmd.append(static_cast<int>(prms.meanR_scaling));
  igpCmd.broadcast();
}

void CLatticeMaster::addPairIG(const CRotThermElasticIGP &prms)
{
  RotIGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.m_kr);
  igpCmd.append(prms.diffusivity);
  igpCmd.broadcast();
}

void CLatticeMaster::addPairIG(const CRotThermFrictionIGP &prms)
{
  RotIGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.k);
  igpCmd.append(prms.mu_s);
  igpCmd.append(prms.mu_d);
  igpCmd.append(prms.k_s);
  igpCmd.append(prms.diffusivity);
  igpCmd.append(prms.dt);
  igpCmd.broadcast();
}

//---
// creation of tagged pair interaction groups
//---


// --- RotFriction ---
class TaggedIGPCommand : public IGPCommand
{
public:
  TaggedIGPCommand(const MpiRankAndComm &globalRankAndComm) : IGPCommand(globalRankAndComm, CMD_ADDTAGPIG)
  {
  }

  void appendIGP(const AIGParam &prms)
  {
    appendTypeAndName(prms);
  }

private:
};

void CLatticeMaster::addTaggedPairIG(
  const CRotFrictionIGP &prms,
  int tag1,
  int mask1,
  int tag2,
  int mask2
)
{
  TaggedIGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.k);
  igpCmd.append(prms.mu_s);
  igpCmd.append(prms.mu_d);
  igpCmd.append(prms.k_s);
  igpCmd.append(prms.dt);
  igpCmd.append(static_cast<int>(prms.scaling));
  igpCmd.append(static_cast<int>(prms.meanR_scaling));
  igpCmd.append(tag1);
  igpCmd.append(mask1);
  igpCmd.append(tag2);
  igpCmd.append(mask2);

  igpCmd.broadcast();
}

// --- NRotFriction ---
void CLatticeMaster::addTaggedPairIG(
  const CFrictionIGP &prms,
  int tag1,
  int mask1,
  int tag2,
  int mask2
)
{
  TaggedIGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.k);
  igpCmd.append(prms.mu);
  igpCmd.append(prms.k_s);
  igpCmd.append(prms.dt);
  igpCmd.append(static_cast<int>(prms.m_scaling));
  igpCmd.append(tag1);
  igpCmd.append(mask1);
  igpCmd.append(tag2);
  igpCmd.append(mask2);
  igpCmd.broadcast();
}

void CLatticeMaster::addTaggedPairIG(
  const CSpringDashpotFrictionIGP &prms,
  int tag1,
  int mask1,
  int tag2,
  int mask2
)
{
  TaggedIGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.youngModulus);
  igpCmd.append(prms.poissonRatio);
  igpCmd.append(prms.cor);
  igpCmd.append(prms.mu);
  igpCmd.append(prms.dt);
  igpCmd.append(tag1);
  igpCmd.append(mask1);
  igpCmd.append(tag2);
  igpCmd.append(mask2);
  igpCmd.broadcast();
}

void CLatticeMaster::addTaggedPairIG(
  const CLinearDashpotIGP &prms,
  int tag1,
  int mask1,
  int tag2,
  int mask2
)
{
  TaggedIGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.m_damp);
  igpCmd.append(prms.m_cutoff);
  igpCmd.append(tag1);
  igpCmd.append(mask1);
  igpCmd.append(tag2);
  igpCmd.append(mask2);

  igpCmd.broadcast();
}


// --- Hertzian Elastic ---
void CLatticeMaster::addTaggedPairIG(
  const CHertzianElasticIGP &prms,
  int tag1,
  int mask1,
  int tag2,
  int mask2
)
{
  TaggedIGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.m_E);
  igpCmd.append(prms.m_nu);
  igpCmd.append(tag1);
  igpCmd.append(mask1);
  igpCmd.append(tag2);
  igpCmd.append(mask2);

  igpCmd.broadcast();
}

// --- Hertzian ViscoElasticFriction ---
void CLatticeMaster::addTaggedPairIG(
  const CHertzianViscoElasticFrictionIGP &prms,
  int tag1,
  int mask1,
  int tag2,
  int mask2
)
{
  TaggedIGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.m_A);
  igpCmd.append(prms.m_E);
  igpCmd.append(prms.m_nu);
  igpCmd.append(prms.mu);
  igpCmd.append(prms.k_s);
  igpCmd.append(prms.dt);
  igpCmd.append(tag1);
  igpCmd.append(mask1);
  igpCmd.append(tag2);
  igpCmd.append(mask2);

  igpCmd.broadcast();
}

// --- Hertzian ViscoElastic ---
void CLatticeMaster::addTaggedPairIG(
  const CHertzianViscoElasticIGP &prms,
  int tag1,
  int mask1,
  int tag2,
  int mask2
)
{
  TaggedIGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.m_A);
  igpCmd.append(prms.m_E);
  igpCmd.append(prms.m_nu);
  igpCmd.append(tag1);
  igpCmd.append(mask1);
  igpCmd.append(tag2);
  igpCmd.append(mask2);

  igpCmd.broadcast();
}

// --- HertzMindlin ---
void CLatticeMaster::addTaggedPairIG(
  const CHertzMindlinIGP &prms,
  int tag1,
  int mask1,
  int tag2,
  int mask2
)
{
  TaggedIGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.m_E);
  igpCmd.append(prms.m_nu);
  igpCmd.append(prms.mu);
  igpCmd.append(prms.dt);
  igpCmd.append(tag1);
  igpCmd.append(mask1);
  igpCmd.append(tag2);
  igpCmd.append(mask2);

  igpCmd.broadcast();
}

// --- HertzMindlinVisco ---
void CLatticeMaster::addTaggedPairIG(
  const CHertzMindlinViscoIGP &prms,
  int tag1,
  int mask1,
  int tag2,
  int mask2
)
{
  TaggedIGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.m_E);
  igpCmd.append(prms.m_nu);
  igpCmd.append(prms.mu);
  igpCmd.append(prms.m_COR);
  igpCmd.append(prms.dt);
  igpCmd.append(tag1);
  igpCmd.append(mask1);
  igpCmd.append(tag2);
  igpCmd.append(mask2);

  igpCmd.broadcast();
}

// --- tagged elastic ---
void CLatticeMaster::addTaggedPairIG(
  const CElasticIGP &prms,
  int tag1,
  int mask1,
  int tag2,
  int mask2
)
{
  TaggedIGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.m_k);
  igpCmd.append(static_cast<int>(prms.m_scaling));
  igpCmd.append(tag1);
  igpCmd.append(mask1);
  igpCmd.append(tag2);
  igpCmd.append(mask2);

  igpCmd.broadcast();
}

void CLatticeMaster::addTaggedPairIG(
  const CRotElasticIGP &prms,
  int tag1,
  int mask1,
  int tag2,
  int mask2
)
{
  TaggedIGPCommand igpCmd(getGlobalRankAndComm());
  igpCmd.appendIGP(prms);
  igpCmd.append(prms.m_kr);
  igpCmd.append(static_cast<int>(prms.m_scaling));
  igpCmd.append(tag1);
  igpCmd.append(mask1);
  igpCmd.append(tag2);
  igpCmd.append(mask2);

  igpCmd.broadcast();
}

/*!
  Remove interaction group. Send name of the interactiongroup to workers

  \param name the name of the interaction group which is to be removed
*/
void CLatticeMaster::removeIG(const std::string& name)
{
  BroadcastCommand cmd(getGlobalRankAndComm(),CMD_REMOVEIG);

  cmd.append(name.c_str());
  cmd.broadcast();
}

//---
// creation of interaction groups between particles and meshes
//---

/*!
 Broadcast command to add a (non-bonded) trimesh IG
*/
class TriMeshIGCommand : public BroadcastCommand
{
public:
  TriMeshIGCommand(const MpiRankAndComm &globalRankAndComm)
    : BroadcastCommand(globalRankAndComm, CMD_ADDTRIMESHIG)
  {
  }
};

/*!
  add(non-bonded) trimesh IG
*/
void CLatticeMaster::addTriMeshIG(const ETriMeshIP &prms)
{
  TriMeshIGCommand cmd(getGlobalRankAndComm());

  cmd.appendTypeAndName(prms);
  cmd.append(prms.getMeshName().c_str());
  cmd.append(prms.k);
  cmd.broadcast();
}

/*!
 Broadcast command to add a (non-bonded) interaction group between particles and
 a 2D mesh
*/
class Mesh2DIGCommand : public BroadcastCommand
{
public:
  Mesh2DIGCommand(const MpiRankAndComm &globalRankAndComm)
    : BroadcastCommand(globalRankAndComm, CMD_ADDMESH2DIG)
  {
  }
};

/*!
  add(non-bonded) interaction group between particles and
 a 2D mesh
*/
void CLatticeMaster::addMesh2DIG(const ETriMeshIP &prms)
{
  Mesh2DIGCommand cmd(getGlobalRankAndComm());

  cmd.appendTypeAndName(prms);
  cmd.append(prms.getMeshName().c_str());
  cmd.append(prms.k);
  cmd.broadcast();
}

class BondedTriMeshIGCommand : public BroadcastCommand
{
public:
  BondedTriMeshIGCommand(const MpiRankAndComm &globalRankAndComm)
    : BroadcastCommand(globalRankAndComm, CMD_ADDBONDEDTRIMESHIG)
  {
  }

  void appendTriMeshPrms(const BTriMeshIP &triMeshPrms)
  {
    append(triMeshPrms.getName().c_str());
    append(triMeshPrms.getMeshName().c_str());
    append(triMeshPrms.k);
    append(triMeshPrms.brk);
  }

  void appendTagBuildPrms(const MeshTagBuildPrms &buildPrms)
  {
    append(buildPrms.getTypeString().c_str());
    append(buildPrms.m_tag);
    append(buildPrms.m_mask);
  }

  void appendGapBuildPrms(const MeshGapBuildPrms &buildPrms)
  {
    append(buildPrms.getTypeString().c_str());
    append(buildPrms.m_maxGap);
  }
};

/*!
  add bonded trimesh IG by particle tag

  \param triMeshPrms the interaction parameters, i.e. k , r_break....
  \param buildPrms the build parameters, i.e. the tag & mask
*/
void CLatticeMaster::addBondedTriMeshIG(const BTriMeshIP &triMeshPrms, const MeshTagBuildPrms &buildPrms)
{
  console.Debug() << "CLatticeMaster::addBondedTriMeshIG with tag build prms\n";
  BondedTriMeshIGCommand bndTriMeshCmd(getGlobalRankAndComm());
  console.Debug() << "CLatticeMaster::addBondedTriMeshIG appending mesh prms to command...\n";
  bndTriMeshCmd.appendTriMeshPrms(triMeshPrms);
  console.Debug() << "CLatticeMaster::addBondedTriMeshIG appending mesh build prms to command...\n";
  bndTriMeshCmd.appendTagBuildPrms(buildPrms);
  console.Debug() << "CLatticeMaster::addBondedTriMeshIG broadcasting...\n";
  bndTriMeshCmd.broadcast();
  console.Debug() << "end CLatticeMaster::addBondedTriMeshIG\n";

}

/*!
  add bonded trimesh IG by distance between particle & mesh

  \param triMeshPrms the interaction parameters, i.e. k , r_break....
  \param buildPrms the build parameters, i.e. the max. dist
*/
void CLatticeMaster::addBondedTriMeshIG(const BTriMeshIP &triMeshPrms, const MeshGapBuildPrms &buildPrms)
{
  console.Debug() << "CLatticeMaster::addBondedTriMeshIG with gap build prms\n";
  BondedTriMeshIGCommand bndTriMeshCmd(getGlobalRankAndComm());
  console.Debug() << "CLatticeMaster::addBondedTriMeshIG appending mesh prms to command...\n";
  bndTriMeshCmd.appendTriMeshPrms(triMeshPrms);
  console.Debug() << "CLatticeMaster::addBondedTriMeshIG appending mesh build prms to command...\n";
  bndTriMeshCmd.appendGapBuildPrms(buildPrms);
  console.Debug() << "CLatticeMaster::addBondedTriMeshIG broadcasting...\n";
  bndTriMeshCmd.broadcast();
  console.Debug() << "end CLatticeMaster::addBondedTriMeshIG\n";
}

/*!
  add bonded interactions with 2d mesh using a gap parameter
*/
void CLatticeMaster::addBondedMesh2DIG(const BMesh2DIP &igp, const MeshGapBuildPrms &buildPrms)
{
  console.Debug() << "CLatticeMaster::addBondedMesh2DIG with Mesh2DGapBuildPrms\n";
  BondedMesh2DIGCommand cmd(getGlobalRankAndComm());
  console.Debug() << "made cmd\n";
  cmd.appendMesh2DParam(igp);
  console.Debug() << "appended IGP\n";
  cmd.appendGapBuildPrms(buildPrms);
  console.Debug() << "appended buildParams\n";
  cmd.broadcast();
  console.Debug() << "end CLatticeMaster::addBondedMesh2DIG\n";
}

/*!
  add bonded interactions with 2d mesh using tag/mask parameters
*/
void CLatticeMaster::addBondedMesh2DIG(const BMesh2DIP &igp, const MeshTagBuildPrms &buildPrms)
{
  console.Debug() << "CLatticeMaster::addBondedMesh2DIG with Mesh2DGapBuildPrms\n";
  BondedMesh2DIGCommand cmd(getGlobalRankAndComm());
  console.Debug() << "made cmd\n";
  cmd.appendMesh2DParam(igp);
  console.Debug() << "appended IGP\n";
  cmd.appendTagBuildPrms(buildPrms);
  console.Debug() << "appended buildParams\n";
  cmd.broadcast();
  console.Debug() << "end CLatticeMaster::addBondedMesh2DIG\n";
}

#if 0
/*!
  add a group of bonded triangle mesh (surface) to particle interactions

    \param name the name of the itneraction group
    \param meshname the name of the mesh
    \param k spring constant
    \param r_break (relative) breaking distance
    \param buildtype the way to initially build the IG
    \param buildparam the parameters for building the IG
*/
void CLatticeMaster::addBondedTriMeshIG(const std::string& name,
                    const std::string& meshname,
                    double k,
                    double r_break,
                    const std::string& buildtype,
                    std::list<OpValue> buildparams)
{
  console.XDebug() << "begin CLatticeMaster::addBondedTriMeshIG\n";
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer pbuffer(m_global_comm);
  CMPIBarrier barrier(m_global_comm);

  // send command
  cmd_buffer.broadcast(CMD_ADDBONDEDTRIMESHIG);

  // send data
  pbuffer.append(name.c_str());
  pbuffer.append(meshname.c_str());
  pbuffer.append(k);
  pbuffer.append(r_break);
  pbuffer.append(buildtype.c_str());
  if(buildtype=="BuildByTag"){
    console.XDebug() << "BuildByTag\n";
    // extract build parameters
    list<OpValue>::iterator iter=buildparams.begin();
    pbuffer.append(int(iter->Content.l)); iter++; // tag
    pbuffer.append(int(iter->Content.l));
    // send params
    pbuffer.broadcast(0);
  } else if (buildtype=="BuildByGap"){
    list<OpValue>::iterator iter=buildparams.begin();
    pbuffer.append(iter->Content.d); // max. gap
    // send params
    pbuffer.broadcast(0);
  } else {
    throw std::runtime_error(std::string("Unknown BTriMeshInteraction build type: ")+buildtype);
  }
  barrier.wait("readBondedTriMeshIG");

  console.XDebug() << "end CLatticeMaster::addBondedTriMeshIG\n";
}

#endif

/*!
  add a triangle mesh

  \param name the name of the mesh
  \param filename the name of the mesh file
*/
void CLatticeMaster::addTriMesh(const string& meshName, const string& fileName)
{
  console.XDebug() << "begin CLatticeMaster::addTriMesh\n";

  readAndDistributeTriMesh(meshName, fileName);

  console.XDebug() << "end CLatticeMaster::addTriMesh\n";
}

/*!
  Read a triangle mesh from a file and distribute the data to the workers.
  The tags on the mesh data are ignored, i.e. the whole file is read as a single mesh.

  \param meshName the name of the mesh
  \param meshFileName the filename
*/
void CLatticeMaster::readAndDistributeTriMesh(const std::string& meshName, const std::string& meshFileName)
{
  TriMeshDataPair meshData = readTriMesh(meshFileName);
  createTriMesh(meshName, meshData.first, meshData.second);
}

/*!
  Read a triangle mesh from a file and distribute the data to the workers.

  \param meshName the name of the mesh
  \param meshFileName the filename
  \param tag the tag in the mesh data determining if a triangle belongs to this mesh or not
*/
void CLatticeMaster::readAndDistributeTriMesh(
     const std::string& meshName,
     const std::string& meshFileName,
     int tag)
{
  TriMeshDataPair meshData = readTriMesh(meshFileName,tag);
  createTriMesh(meshName, meshData.first, meshData.second);
}

/*!
  Read triangle mesh data from file and return data as <vector of nodes,vector of triangles> pair.
  Tags are ignored, i.e. the whole file is read as a single mesh.

  \param meshfilename the filename
*/
CLatticeMaster::TriMeshDataPair CLatticeMaster::readTriMesh(const std::string& meshfilename)
{
  console.XDebug() << "CLatticeMaster::readTriMesh(" << meshfilename << ")\n";
  // buffers
  MeshNodeDataVector node_send_buffer;
  MeshTriDataVector tri_send_buffer;

  // mesh reader
  MeshReader reader(meshfilename);

  // --- Nodes ---
  MeshReader::NodeIterator &niter=reader.getNodeIterator();
  // read nodes into vector
  while(niter.hasNext()){
    node_send_buffer.push_back(niter.next());
  }

  // --- Triangles ---
  MeshReader::TriIterator &titer=reader.getTriIterator();
  // read triangles into vector
  while(titer.hasNext()){
    tri_send_buffer.push_back(titer.next());
    // debug
    // MeshTriData last=tri_send_buffer.back();
    // console.XDebug() << "got triangle data - id: " << last.id << " tag: " << last.tag << "\n";
  }
  return TriMeshDataPair(node_send_buffer, tri_send_buffer);
}

/*!
  read a triangle mesh from a file and distribute the data to the workers

  \param meshfilename the filename
  \param tag the tag in the mesh data determining if a triangle belongs to this mesh or not
*/
CLatticeMaster::TriMeshDataPair CLatticeMaster::readTriMesh(const std::string& meshfilename,int tag)
{
  console.XDebug() << "CLatticeMaster::readTriMesh(" << meshfilename << ")\n";
  // buffers
  vector<MeshNodeData> node_buffer;
  vector<MeshTriData> tri_buffer;

  // mesh reader
  MeshReader reader(meshfilename);


  // --- Nodes ---
  MeshReader::NodeIterator &niter=reader.getNodeIterator();
  // read nodes into vector
  while(niter.hasNext()){
    node_buffer.push_back(niter.next());
  }

  // --- Triangles ---
  MeshReader::TriIterator &titer=reader.getTriIterator();
  // read triangles into vector
  while(titer.hasNext()){
    MeshTriData mtd=titer.next();
    if(mtd.tag==tag){
      tri_buffer.push_back(mtd);
    }
  }

  console.XDebug() << "end CLatticeMaster::readTriMesh\n";
  return TriMeshDataPair(node_buffer, tri_buffer);
}

/*!
  create a triangle mesh from a vector of nodes and a vector of triangles

  \param meshName the name assigned to the triangle mesh
  \param node_send_buffer the vector of nodes
  \param tri_send_buffer the vector of triangles
*/
void CLatticeMaster::createTriMesh(
  const std::string &meshName,
  const MeshNodeDataVector &node_send_buffer,
  const MeshTriDataVector &tri_send_buffer
)
{
  console.XDebug() << "CLatticeMaster::distributeTriMesh: enter\n";
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer pbuffer(m_global_comm);

  // send command
  cmd_buffer.broadcast(CMD_ADDTRIMESH);

  pbuffer.append(meshName.c_str());
  pbuffer.broadcast(0);

  CMPIBarrier barrier(m_global_comm);
  // distribute the node vector
  m_tml_global_comm.broadcast_cont_packed(node_send_buffer);
  // distribute the triangle vector
  m_tml_global_comm.broadcast_cont_packed(tri_send_buffer);

  barrier.wait("distributeTriMesh");
  console.XDebug() << "CLatticeMaster::distributeTriMesh: exit\n";
}


/*!
  add a 2D mesh

  \param name the name of the mesh
  \param filename the name of the mesh file
  \param tag the tag of the edges that are included into the mesh
*/
void CLatticeMaster::addMesh2D(const string& name, const string& filename, int tag)
{
  console.XDebug() << "CLatticeMaster::addMesh2D( " << name << " , " << filename << ")\n";

  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  CVarMPIBuffer pbuffer(m_global_comm);

  // send command
  cmd_buffer.broadcast(CMD_ADDMESH2D);

  pbuffer.append(name.c_str());
  pbuffer.broadcast(0);
  readAndDistributeMesh2D(filename,tag);

  console.XDebug() << "end CLatticeMaster::addMesh2D\n";
}

/*!
  read a 2D mesh from a file and distribute the data to the workers

  \param meshfilename the filename
  \param tag the tag of the edges that are included into the mesh
*/
void CLatticeMaster::readAndDistributeMesh2D(const std::string& meshfilename,int tag)
{
  console.XDebug() << "CLatticeMaster::readAndDistributeMesh2D(" << meshfilename << ")\n";
  CMPIBarrier barrier(m_global_comm);
  vector<MeshNodeData2D> node_send_buffer;
  vector<MeshEdgeData2D> edge_send_buffer;

  // mesh reader
  Mesh2DReader reader(meshfilename);

  // --- Nodes ---
  Mesh2DReader::NodeIterator &niter=reader.getNodeIterator();
  console.XDebug() << "Node reader initialized\n";
  // read nodes into vector
  while(niter.hasNext()){
    node_send_buffer.push_back(niter.next());
  }

  // --- Edges ---
  Mesh2DReader::EdgeIterator &eiter=reader.getEdgeIterator();
  console.XDebug() << "Edge reader initialized\n";
  // read triangles into vector
  while(eiter.hasNext()){
    MeshEdgeData2D med=eiter.next();
    if(med.tag==tag){
      edge_send_buffer.push_back(med);
    }
  }

  // distribute the node vector
  m_tml_global_comm.broadcast_cont_packed(node_send_buffer);
  // distribute the triangle vector
  m_tml_global_comm.broadcast_cont_packed(edge_send_buffer);

  barrier.wait("readAndDistributeMesh2D");
  console.XDebug() << "end CLatticeMaster::readAndDistributeMesh2D\n";
}

class SIGCommand : public BroadcastCommand
{
public:
  SIGCommand(const MpiRankAndComm &globalRankAndComm)
    : BroadcastCommand(globalRankAndComm, CMD_ADDSIG)
  {
  }

  void appendGravityIGP(const esys::lsm::GravityIGP &gravityIGP)
  {
    append(gravityIGP.getTypeString().c_str());
    packInto(gravityIGP);
  }

  void appendBuoyancyIGP(const esys::lsm::BuoyancyIGP &buoyancyIGP)
  {
    append(buoyancyIGP.getTypeString().c_str());
    packInto(buoyancyIGP);
  }
};

void CLatticeMaster::addSingleIG(const esys::lsm::GravityIGP &gravityIGP)
{
  SIGCommand sigCmd(getGlobalRankAndComm());
  sigCmd.appendGravityIGP(gravityIGP);
  sigCmd.broadcast();
}

void CLatticeMaster::addSingleIG(const esys::lsm::BuoyancyIGP &buoyancyIGP)
{
  SIGCommand sigCmd(getGlobalRankAndComm());
  sigCmd.appendBuoyancyIGP(buoyancyIGP);
  sigCmd.broadcast();
}

#if 0
/*!
  add body force

  \param type the type of force, i.e. Gravity,...
  \param prms the interaction parameters
*/
void CLatticeMaster::addBodyForce(const string& type, const esys::lsm::BodyForceIGP &prms)
{
  console.XDebug() << "begin CLatticeMaster::addBodyForce\n";
  if (type == "Gravity") {
    console.XDebug()<<"Body force type is " << type << "\n";

    CVarMPIBuffer pbuffer(MPI_COMM_WORLD,20);
    CMPIBarrier barrier(MPI_COMM_WORLD);
    CMPILCmdBuffer cmd_buffer(MPI_COMM_WORLD,m_global_rank);

    //send the command
    cmd_buffer.broadcast(CMD_ADDSIG);
    // console.Debug()<<"adding Gravity IG\n";
    pbuffer.append(type.c_str());
    prms.packInto(&pbuffer);
    pbuffer.broadcast(0);

    barrier.wait("addBodyForce()");
  }
  else {
    throw std::runtime_error(
      std::string("Unknown body force type: ")
      +
      type
    );
  }
  console.XDebug() << "end CLatticeMaster::addBodyForce\n";
}

void CLatticeMaster::addSingleIG(const string& type, const string& name, list<OpValue> params)
{
  console.XDebug()<<"CLatticeMaster::addSingleIG\n";

  if(type=="Gravity"){

    list<OpValue>::iterator iter=params.begin();
    double x = iter->Content.d;
    iter++;
    double y = iter->Content.d;
    iter++;
    double z = iter->Content.d;

    addBodyForce(type, esys::lsm::BodyForceIGP(name, Vec3(x, y, z)));
  }else{
    throw std::runtime_error(
      std::string("Unknown body force type: ")
      +
      type
    );
  }
  console.XDebug()<<"end CLatticeMaster::addSingleIG\n";
}
#endif

class DampingCommand : public BroadcastCommand
{
public:
  DampingCommand(const MpiRankAndComm &globalRankAndComm)
    : BroadcastCommand(globalRankAndComm, CMD_ADDDAMP)
  {
  }
};

void CLatticeMaster::addDamping(const CDampingIGP &dampingIGP)
{
  DampingCommand dampCmd(getGlobalRankAndComm());
  dampCmd.append(dampingIGP.getTypeString().c_str());
  dampCmd.packInto(dampingIGP);
  dampCmd.broadcast();
}

void CLatticeMaster::addDamping(const CLocalDampingIGP &dampingIGP)
{
  DampingCommand dampCmd(getGlobalRankAndComm());
  dampCmd.append(dampingIGP.getTypeString().c_str());
  dampCmd.packInto(dampingIGP);
  dampCmd.broadcast();
}

void CLatticeMaster::addDamping(const ABCDampingIGP &dampingIGP)
{
  DampingCommand dampCmd(getGlobalRankAndComm());
  dampCmd.append(dampingIGP.getTypeString().c_str());
  dampCmd.packInto(dampingIGP);
  dampCmd.broadcast();
}

#if 0
/*!
  Add damping to a lattice. Send parameters to slaves
*/
void CLatticeMaster::addDamping(const string& type,list<OpValue> params)
{
  console.XDebug()<<"CLatticeMaster::addDamping\n";

  CVarMPIBuffer pbuffer(m_global_comm,20);
  CMPIBarrier barrier(m_global_comm);
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);

  //send the command
  cmd_buffer.broadcast(CMD_ADDDAMP);
  if(type=="Damping"){
    // console.Debug()<<"adding Damping\n";
    pbuffer.append("Damping");
    CDampingIGP igp;
    list<OpValue>::iterator iter=params.begin();
    igp.setName("damp");
    igp.setVRef(Vec3(0.0,0.0,0.0));
    igp.setVisc(iter->Content.d);
    console.XDebug()<< "visc : " << iter->Content.d;
    iter++;
    igp.setTimeStep(iter->Content.d);
    console.XDebug()<< "dt : " << iter->Content.d;
    iter++;
    igp.setMaxIter(iter->Content.l);
    igp.packInto(&pbuffer);
    pbuffer.broadcast(0);
  }else{
    console.Error()<<"trying to add Damping of unknown type\n";
  }
  barrier.wait("addDamping()");
  console.XDebug()<<"end CLatticeMaster::addDamping\n";
}

#endif

/*!
  set IG s1 as excluding IG in s2, i.e. an two particles interacting in s1
  can't interact in s2
*/
void CLatticeMaster::addExIG(const string& s1,const string& s2)
{
  console.XDebug()<<"CLatticeMaster::addExIG\n";

  CVarMPIBuffer pbuffer(m_global_comm,20);
  CMPIBarrier barrier(m_global_comm);
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);

  //send the command
  cmd_buffer.broadcast(CMD_EXIG);
  pbuffer.append(s1.c_str());
  pbuffer.append(s2.c_str());
  pbuffer.broadcast(0);
  barrier.wait("addExIG()");
  console.XDebug()<<"end CLatticeMaster::addExIG\n";
  // update interactions - relevant in case the exclusion is added _during_ a run
  updateInteractions();
}

#if 0
/*!
  Get the list of face ids for a given mesh, i.e. edge ids in 2D, triangle ids in 3D

  \param meshname
*/
void CLatticeMaster::getMeshFaceReferences(const string& meshname)
{
  console.XDebug()<<"CLatticeMaster::getMeshFaceReferences( " << meshname << ")\n";
  GetFaceRefCommand cmd(getGlobalRankAndComm(),meshname);
  // broadcast command
  cmd.broadcast();

  // receive data (multimap)
  multimap<int,int> ref_mmap;
  m_tml_global_comm.gather(ref_mmap);
  // collate into set
  set<int> ref_set; //== this is the set of node ids ==
  for(multimap<int,int>::iterator iter=ref_mmap.begin();
      iter!=ref_mmap.end();
      iter++){
    ref_set.insert(iter->second);
  }
  console.XDebug()<<"end CLatticeMaster::getMeshFaceReferences()\n";
}

/*!
  Get the list of stresses for a given mesh.

  \param meshname
*/
void CLatticeMaster::getMesh2DEdgeStress(const string& meshname)
{
  console.XDebug()<<"CLatticeMaster::getMesh2DEdgeStress( " << meshname << ")\n";
  multimap<int,pair<int,Vec3> > temp_mm;
  map<int,Vec3> m_data; //=== map of id, value ===

  BroadcastCommand cmd(getGlobalRankAndComm(), CMD_GETMESH2DSTRESS);
  cmd.append(meshname.c_str());
  cmd.broadcastCommand();


  // get data from slaves
  m_tml_global_comm.gather(temp_mm);

  // add data together
  for(multimap<int,pair<int,Vec3> >::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    if(m_data.find((iter->second).first)==m_data.end()){ // id not in m_data -> insert
      m_data.insert(iter->second);
    } else { // id is in m_data -> add
      m_data[(iter->second).first]+=(iter->second).second;
    }
  }

  console.XDebug()<<"end CLatticeMaster::getMesh2DEdgeStress()\n";
}
#endif

/*!
  Set local console verbosity and send command to workers to set console verbosity there
*/
void CLatticeMaster::setVerbosity(int verbose)
{
  console.SetVerbose(verbose);

  // workers
  BroadcastCommand cmd(getGlobalRankAndComm(), CMD_SETVERBOSITY);
  cmd.append(verbose);
  cmd.broadcast();

}

/*!
  Set local console filename and send command to workers to set console parameters there
*/
void CLatticeMaster::setConsoleFilename(const string& fname)
{
    console.SetFilename(fname);

    // workers
    BroadcastCommand cmd(getGlobalRankAndComm(), CMD_SETCONSOLEFNAME);
    cmd.append(fname.c_str());
    cmd.broadcast();
}

/*!
  Set local console buffer size & buffering mode and send command to
  workers to set console parameters there
*/
void CLatticeMaster::setConsoleBuffered(unsigned int bsize)
{
    console.SetBuffered(bsize);

    // workers
    BroadcastCommand cmd(getGlobalRankAndComm(), CMD_SETCONSOLEBUFF);
    cmd.append(int(bsize));
    cmd.broadcast();
}

/*!
  initialize local console and send command to workers initialize consoles there

  \param filename the base name of the output file
  \param bufflen the length of the internal buffer in the console 0-> no buffering
*/
void CLatticeMaster::initializeConsole(const string& filename, int bufflen)
{
  // workers
  BroadcastCommand cmd(getGlobalRankAndComm(), CMD_INITCONSOLE);
  cmd.append(filename.c_str());
  cmd.append(bufflen);
  cmd.broadcast();

  //master
  console.Initialize(filename);
  console.SetBuffered(bufflen);
}

/*!
  Set a named parameter of a all interactions in a given interaction group to
  a specific value.
   
  \param igname the name of the affected interaction group
  \param pname the name of the parameter
  \param val the value 
*/
void CLatticeMaster::setInteractionParameter(const string& igname,const string& pname,double val)
{
    console.XDebug()<<"CLatticeMaster::setInteractionParameter( " << igname << " , " << pname << " , " << val << ")\n";
    
    // pack parameters and send to workers
    BroadcastCommand cmd(getGlobalRankAndComm(), CMD_SETINTERACTIONPARAMS);
    cmd.append(igname.c_str());
    cmd.append(pname.c_str());
    cmd.append(val);
    cmd.broadcast();
    
    console.XDebug()<<"end CLatticeMaster::setInteractionParameter\n";
}
