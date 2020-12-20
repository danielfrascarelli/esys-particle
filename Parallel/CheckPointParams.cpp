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


#include "Parallel/CheckPointParams.h"
#include "Parallel/MpiInfo.h"
#include "Parallel/mpibuf.h"

#include <sstream>
#include <memory>

CheckPointParams::CheckPointParams() : m_fileNamePrefix(), m_time(-1), m_rank(-1), m_prec(12)
{
}

CheckPointParams::CheckPointParams
(
  const std::string &fileNamePrefix,
  int time,
  int rank,
  int prec
) : m_fileNamePrefix(fileNamePrefix),
    m_time(time),
    m_rank(rank),
    m_prec(prec)
{}

CheckPointParams::~CheckPointParams()
{
}

std::string CheckPointParams::getFileName() const
{
  std::stringstream name;
  name << m_fileNamePrefix << "_t=" << m_time << "_" << m_rank << ".txt";
  
  return name.str();
}

void CheckPointParams::packInto(AMPIBuffer *pMpiBuff) const
{
  pMpiBuff->append(m_fileNamePrefix.c_str());
  pMpiBuff->append(m_time);
  pMpiBuff->append(m_prec);
}

CheckPointParams CheckPointParams::unpackFrom(AMPIBuffer *pMpiBuff, MPI_Comm mpiComm)
{
  CheckPointParams prms;

  prms.m_fileNamePrefix = pMpiBuff->pop_string();
  prms.m_time = pMpiBuff->pop_int();
  prms.m_prec = pMpiBuff->pop_int();
  prms.m_rank = MpiInfo(mpiComm).rank();

  return prms;
}
