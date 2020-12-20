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
#include "Model/ShortBondedInteractionCpData.h"
#include "Model/ShortBondedInteraction.h"

// --- STL includes ---
#include <utility>

using std::pair;

/*!
  default constructor
*/
ShortBondedInteractionCpData::ShortBondedInteractionCpData()
{}

/*!
  constructor

  \param p1id id of the 1st particle
  \param p2id id of the 2nd particle
  \param itag interaction tag
  \param r0 equilibrium distance
*/
ShortBondedInteractionCpData::ShortBondedInteractionCpData(int p1id,int p2id,int itag,double r0)
  : BondedInteractionCpData(p1id,p2id,itag)
{
  m_r0=r0;
}

/*!
  construct directly from CShortBondedInteraction

  \param SBI the interaction
*/
ShortBondedInteractionCpData::ShortBondedInteractionCpData(const CShortBondedInteraction &SBI)
{
  pair<int,int> pids=SBI.getPairID();
  int itag=SBI.getTag();
  double r0=SBI.getEquiDist();
  set(pids.first,pids.second,itag); // base 
  m_r0=r0;
}
 
/*!
  write data to output stream

  \param ost the output stream
*/
void ShortBondedInteractionCpData::saveCheckPointData(ostream& ost)
{
  const char delim = ' ';
  BondedInteractionCpData::saveCheckPointData(ost);
  ost << delim << m_r0;
}

/*!
  read data from input stream

  \warning not implemented
*/
void ShortBondedInteractionCpData::loadCheckPointData(istream&)
{}
