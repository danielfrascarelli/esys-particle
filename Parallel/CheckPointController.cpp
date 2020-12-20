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


#include "Parallel/CheckPointController.h"
#include "Parallel/mpicmdbuf.h"
#include "Parallel/mpibarrier.h"
#include "Parallel/sublattice_cmd.h"
#include "Parallel/MpiInfo.h"
#include "Parallel/CheckPointParams.h"
#include "Parallel/mpivbuf.h"
#include "Parallel/CheckPointInfo.h"

//--- TML includes ---
#include "tml/comm/comm.h"

#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>
#include <sstream>

using std::multimap;
using std::vector;
using std::string;
using std::istream_iterator;
using std::back_inserter;
using std::ostringstream;

/*!
  default constructor
*/
CheckPointController::CheckPointController()
  : m_mpiComm(MPI_COMM_WORLD),
    m_fileNamePrefix(""),
    m_beginTime(-1),
    m_endTime(-1),
    m_timeInterval(-1),
    m_geoInfo(),
    m_numTimeSteps(0),
    m_timeStepSize(0.0),
    m_writeThroughMaster(false),
    m_precision(12)
{
  m_spatialDomainHasBeenSet=false;
}

CheckPointController::CheckPointController(
  const std::string &fileNamePrefix,
  int beginTime,
  int endTime,
  int timeInterval,
  bool writeThroughMaster
) : m_mpiComm(MPI_COMM_WORLD),
    m_fileNamePrefix(fileNamePrefix),
    m_beginTime(beginTime),
    m_endTime(endTime),
    m_timeInterval(timeInterval),
    m_geoInfo(),
    m_numTimeSteps(0),
    m_timeStepSize(0.0),
    m_writeThroughMaster(writeThroughMaster),
    m_precision(12)
{
  m_spatialDomainHasBeenSet=false;
}

CheckPointController::~CheckPointController()
{
}

/*!
  Issue checkpointing command, i.e. broadcast command via MPI to
  worker processes, if currentTime is a time at which a checkpoint 
  needs to be taken

  \param currentTime the current time step
*/
void CheckPointController::performCheckPoint(int currentTime)
{
  if (isCheckPoint(currentTime))
  {
    if(m_writeThroughMaster){
      issueCheckPointCmdWTM(currentTime);
    } else {
      issueCheckPointCmd(currentTime);
    }
  }
}

/*!
  Issue snapshot command, i.e. broadcast command via MPI to
  worker processes, if currentTime is a time at which a snapshot 
  needs to be taken

  \param currentTime the current time step
*/
void CheckPointController::performSnapShot(int currentTime)
{
  if (isCheckPoint(currentTime))
  {
    issueSnapShotCmd(currentTime);
  }
}

MPI_Comm CheckPointController::getMpiComm() const
{
  return m_mpiComm;
}

void CheckPointController::setMpiComm(MPI_Comm mpiComm)
{
  m_mpiComm = mpiComm;
}

/*!
  Broadcast checkpointing command to workers and write info file (*_0.txt)
  
  \param currentTime the current time step
*/
void CheckPointController::issueCheckPointCmd(int currentTime)
{
  CMPILCmdBuffer cmd_buffer(getMpiComm(), MpiInfo(getMpiComm()).rank());
  CMPIBarrier barrier(getMpiComm());

  cmd_buffer.broadcast(CMD_SAVECHECKPOINT);
  
  CVarMPIBuffer buffer(getMpiComm());

  // send check-point parameters
  CheckPointParams checkPointParams(m_fileNamePrefix, currentTime, MpiInfo(getMpiComm()).rank(),m_precision);
  checkPointParams.packInto(&buffer);
  buffer.broadcast(0);
  buffer.clear();

  std::ofstream oStream(checkPointParams.getFileName().c_str());
  esys::lsm::CheckPointInfo cpInfo;

  cpInfo.setLatticeDataFiles(getLatticeDataFiles(currentTime));
  cpInfo.setGeometryInfo(getGeometryInfo());
  cpInfo.setTimeStep(currentTime);
  cpInfo.setTimeStepSize(getTimeStepSize());
  cpInfo.setNumTimeSteps(getNumTimeSteps());
  cpInfo.write(oStream);

  oStream.close();

  barrier.wait("CheckPoint");
}

/*!
  checkpointing with writing through master
  
  \param currentTime the current time step
*/
void CheckPointController::issueCheckPointCmdWTM(int currentTime)
{
  CMPILCmdBuffer cmd_buffer(getMpiComm(), MpiInfo(getMpiComm()).rank());
  CMPIBarrier barrier(getMpiComm());

  cmd_buffer.broadcast(CMD_SAVECHECKPOINTWTM);
  
  CVarMPIBuffer buffer(getMpiComm());

  // send check-point parameters
  CheckPointParams checkPointParams(m_fileNamePrefix, currentTime, MpiInfo(getMpiComm()).rank(),m_precision);
  checkPointParams.packInto(&buffer);
  buffer.broadcast(0);
  buffer.clear();

  std::ofstream oStream(checkPointParams.getFileName().c_str());
  esys::lsm::CheckPointInfo cpInfo;

  cpInfo.setLatticeDataFiles(getLatticeDataFiles(currentTime));
  cpInfo.setGeometryInfo(getGeometryInfo());
  cpInfo.setTimeStep(currentTime);
  cpInfo.setTimeStepSize(getTimeStepSize());
  cpInfo.setNumTimeSteps(getNumTimeSteps());
  cpInfo.write(oStream);

  oStream.close();

  barrier.wait("CheckPoint_1");
  
  // receive data from workers
  multimap<int,char> worker_data;
  TML_Comm global_comm(MPI_COMM_WORLD);
    
  global_comm.gather(worker_data);
  
  console.Debug() << "worker data size: " << worker_data.size() << "\n";

  barrier.wait("CheckPoint_2");
}

/*!
  save meta-data and issue snapshot command to worker processes
  N.B.: the "binary" flag is forced to false of non-restart snapshots

  \param currentTime The current time step. 
*/
void CheckPointController::issueSnapShotCmd(int currentTime)
{
  CMPILCmdBuffer cmd_buffer(getMpiComm(), MpiInfo(getMpiComm()).rank());
  CMPIBarrier barrier(getMpiComm());

  cmd_buffer.broadcast(CMD_SAVESNAPSHOT);
  
  CVarMPIBuffer buffer(getMpiComm());

  // send check-point parameters
  CheckPointParams checkPointParams(m_fileNamePrefix, currentTime, MpiInfo(getMpiComm()).rank(),false);
  checkPointParams.packInto(&buffer);
  buffer.broadcast(0);
  buffer.clear();

  ofstream oStream(checkPointParams.getFileName().c_str());
 
  esys::lsm::CheckPointInfo cpInfo;
  cpInfo.setLatticeDataFiles(getLatticeDataFiles(currentTime));
  cpInfo.setGeometryInfo(getGeometryInfo());
  cpInfo.setTimeStep(currentTime);
  cpInfo.setTimeStepSize(getTimeStepSize());
  cpInfo.setNumTimeSteps(getNumTimeSteps());
  cpInfo.write(oStream);

  oStream.close();

  barrier.wait("SnapShot");
}

/*!
  read meta-data and issue checkpoint loading command to worker processes

  \param metafile_name the name of the file with the meta-data
*/
void CheckPointController::issueCheckPointLoadingCmd(const std::string& metafile_name)
{
  console.Debug() << "CheckPointController::issueCheckPointLoadingCmd()\n";

  CMPILCmdBuffer cmd_buffer(getMpiComm(), MpiInfo(getMpiComm()).rank());
  CMPIBarrier barrier(getMpiComm());

  // send command to workers
  cmd_buffer.broadcast(CMD_LOADCHECKPOINT);
  
  CVarMPIBuffer buffer(getMpiComm());

  // read meta-data
  std::ifstream iStream(metafile_name.c_str());
  float dummy,xmax,ymax,zmax,xmin,ymin,zmin,geo_version,version;
  vector<string> filenames;
  string dummystring;
  int currentTime;

  iStream >> dummystring >> version;
  iStream >> dummy >> dummy >> currentTime;
  iStream >> dummystring >> geo_version;

  // get bounding box
  iStream >> dummystring;
  iStream >> xmin >> ymin >> zmin >> xmax >> ymax >> zmax ;

  // ignore periodic bdry
  iStream >> dummystring >> dummy >> dummy >> dummy;
  

  // dimension info (2D/3D)
  iStream >> dummystring >> dummystring;
  

  // get file names
  copy(istream_iterator<string>(iStream),istream_iterator<string>(),back_inserter(filenames));
  
  iStream.close();

  // --- generate checkpoint parameters
  // get nr. of filenames 
  int nr_fnames=filenames.size();
  console.Debug() << "read " <<  nr_fnames << " file names\n";

  // extract prefix from filenames
  string::size_type p=(filenames.begin())->find("_t");
  string prefix=(filenames.begin())->substr(0,p);
  console.Debug() << "found file prefix: " << prefix  << " time step: " << currentTime <<"\n";
  // construct check-point parameters and distribute to workers 
  CheckPointParams checkPointParams(prefix, currentTime, MpiInfo(getMpiComm()).rank(),m_precision);
  checkPointParams.packInto(&buffer);
  buffer.broadcast(0);
  buffer.clear();

  barrier.wait("loadCheckPoint");
}



/*!
  Generate the filename for the checkpoint file to be written by a given worker 

  \param fileNamePrefix global filename prefix 
  \param timeStep current time step
  \param rank the MPI rank of the worker
*/
std::string CheckPointController::getLatticeDataFileName(const std::string &fileNamePrefix, int timeStep, int rank,bool bin)
{
  return CheckPointParams(fileNamePrefix, timeStep, rank,bin).getFileName();
}


/*!
  Generate the filenames for the checkpoint files to be written by all  workers 

  \param timeStep current time step
  \param size number of workers
*/
esys::lsm::StringVector CheckPointController::getLatticeDataFiles(int timeStep, int size)
{
  esys::lsm::StringVector fileNames;
  for (int rank = 1; rank < size; rank++)
  {
    fileNames.push_back(getLatticeDataFileName(m_fileNamePrefix, timeStep, rank));
  }
  return fileNames;
}

/*!
  Generate the filenames for the checkpoint files to be written by all  workers 

  \param timeStep current time step
*/
esys::lsm::StringVector CheckPointController::getLatticeDataFiles(int timeStep)
{
  MpiInfo mpiInfo;
  
  return getLatticeDataFiles(timeStep, MpiInfo(getMpiComm()).size());
}

bool CheckPointController::isCheckPoint(int time)
{
  return
    (
      (
        (time >= m_beginTime)
        &&
        (time <= m_endTime)
        &&
        (((time - m_beginTime) % m_timeInterval) == 0)
      )
      ||
      (time == m_endTime)
    );
}

void CheckPointController::setCheckPointParams(
  const std::string &fileNamePrefix,
  int beginTime,
  int endTime,
  int timeInterval,
  bool writeThroughMaster,
  int precision
)
{
  m_fileNamePrefix = fileNamePrefix;
  m_beginTime      = beginTime;
  m_endTime        = endTime;
  m_timeInterval   = timeInterval;
  m_writeThroughMaster = writeThroughMaster;
  m_precision = precision;
}

void CheckPointController::set_is2d(bool do2d)
{
  m_geoInfo.set_is2d(do2d);
}

void CheckPointController::setLsmGeoVersion(float version)
{
  m_geoInfo.setLsmGeoVersion(version);
}

void CheckPointController::setPeriodicDimensions(esys::lsm::BoolVector periodicDimensions)
{
  m_geoInfo.setPeriodicDimensions(periodicDimensions);
}

void CheckPointController::setSpatialDomain(const esys::lsm::BoundingBox &bbox)
{
  m_geoInfo.setBBox(bbox.getMinPt(),bbox.getMaxPt());
  m_spatialDomainHasBeenSet=true;
}

void CheckPointController::setGeometryInfo(const esys::lsm::GeometryInfo &geoInfo)
{
  m_geoInfo = geoInfo;
}

esys::lsm::GeometryInfo CheckPointController::getGeometryInfo() const
{
  return m_geoInfo;
}

int CheckPointController::getNumTimeSteps() const
{
  return m_numTimeSteps;
}

void CheckPointController::setNumTimeSteps(int numTimeSteps)
{
  m_numTimeSteps = numTimeSteps;
}

double CheckPointController::getTimeStepSize() const
{
  return m_timeStepSize;
}

void CheckPointController::setTimeStepSize(double timeStepSize)
{
  m_timeStepSize = timeStepSize;
}

/*!
  return true if the spatial domain has been set, false otherwise
*/
bool CheckPointController::spatialDomainHasBeenSet() const
{
  return m_spatialDomainHasBeenSet;
}
