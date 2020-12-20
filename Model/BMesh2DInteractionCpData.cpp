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

#include "Model/BMesh2DInteractionCpData.h"
#include "Model/BMesh2DInteraction.h"

BMesh2DInteractionCpData::BMesh2DInteractionCpData() : m_ap(), m_pid(-1), m_tid(-1)
{}

BMesh2DInteractionCpData::BMesh2DInteractionCpData(const BEdge2DInteraction& bmi)
{
  m_pid=bmi.getPid();
  m_tid=bmi.getTid();
  m_ap=bmi.getAP();
}

BMesh2DInteractionCpData::BMesh2DInteractionCpData(int pid,int tid)
{
  m_pid=pid;
  m_tid=tid;
}

void BMesh2DInteractionCpData::set(const BEdge2DInteraction& bmi)
{
  m_pid=bmi.getPid();
  m_tid=bmi.getTid();
  m_ap=bmi.getAP();

}

void BMesh2DInteractionCpData::set(int pid,int tid)
{
  m_pid=pid;
  m_tid=tid;
}

int BMesh2DInteractionCpData::getPID()
{
  return m_pid;
}

int BMesh2DInteractionCpData::getTID()
{
  return m_tid;
}
 
void BMesh2DInteractionCpData::saveSnapShotData(ostream& ost)
{
  const char delim = ' ';
  ost << m_pid << delim << m_tid << delim << m_ap ;
}

void BMesh2DInteractionCpData::saveCheckPointData(ostream& ost)
{
  const char delim = ' ';
  ost << m_pid << delim << m_tid ;
}

void BMesh2DInteractionCpData::loadCheckPointData(istream &ist)
{
  ist >> m_pid >> m_tid;
}
