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


#include "Parallel/CheckPointer.h"
#include "Parallel/CheckPointParams.h"
#include "Parallel/CheckPointable.h"
#include "Parallel/mpivbuf.h"
#include "Parallel/mpibarrier.h"
#include "Foundation/console.h"

//--- TML includes ---
#include "tml/comm/comm.h"

#include <fstream>
#include <sstream>

using esys::lsm::CheckPointable;

CheckPointer::CheckPointer(CheckPointable &checkPointable, MPI_Comm mpiComm)
  : m_pCheckPointable(&checkPointable),
    m_mpiComm(mpiComm)
{
}

CheckPointer::~CheckPointer()
{
}

MPI_Comm CheckPointer::getMpiComm() const
{
  return m_mpiComm;
}

void CheckPointer::setMpiComm(MPI_Comm mpiComm)
{
  m_mpiComm = mpiComm;
}

void CheckPointer::saveRestartable()
{
  CVarMPIBuffer buffer(getMpiComm());

  // get check-point parameters
  buffer.receiveBroadcast(0);
  CheckPointParams checkPointParams = CheckPointParams::unpackFrom(&buffer, getMpiComm());

  // Write the check-point info to file.  
  std::ofstream oStream(checkPointParams.getFileName().c_str());
  // set output precision
  console.Debug() << "output precision: " << checkPointParams.getPrecision() << "\n";

  oStream.precision(checkPointParams.getPrecision());

  m_pCheckPointable->saveCheckPointData(oStream);
  oStream.close();
}

void CheckPointer::saveDump()
{
  CVarMPIBuffer buffer(getMpiComm());

  // get check-point parameters
  buffer.receiveBroadcast(0);
  CheckPointParams checkPointParams = CheckPointParams::unpackFrom(&buffer, getMpiComm());

  // Write the check-point info to file.  
  std::ofstream oStream(checkPointParams.getFileName().c_str());
  m_pCheckPointable->saveSnapShotData(oStream);
  oStream.close();
}

void CheckPointer::saveThroughMaster(TML_Comm& comm)
{
  CVarMPIBuffer buffer(getMpiComm());
  CMPIBarrier barrier(getMpiComm());  

  console.Debug() << "CheckPointer::saveThroughMaster" << "\n";
  // get check-point parameters
  buffer.receiveBroadcast(0);
  CheckPointParams checkPointParams = CheckPointParams::unpackFrom(&buffer, getMpiComm());

  // Write the check-point info to string stream.  
  std::ostringstream oStream;
  save(oStream);
  console.Debug() << "string length : " << oStream.str().size() << "\n";
  string str_data=oStream.str();

  
  barrier.wait("CheckPoint_1");
  comm.send_gather(str_data,0);
}

void CheckPointer::loadCheckPoint()
{
  CVarMPIBuffer buffer(getMpiComm());

  // get check-point parameters
  buffer.receiveBroadcast(0);
  CheckPointParams checkPointParams = CheckPointParams::unpackFrom(&buffer, getMpiComm());

  // Write the check-point info to file.  
  std::ifstream iStream(checkPointParams.getFileName().c_str());
  m_pCheckPointable->loadCheckPointData(iStream);
  iStream.close();
}



void CheckPointer::save(std::ostream &oStream)
{
  m_pCheckPointable->saveCheckPointData(oStream);
}
