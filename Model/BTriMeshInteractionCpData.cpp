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

#include "Model/BTriMeshInteractionCpData.h"
#include "Model/BTriangleInteraction.h"

/*!
  Default constructor. Particle and Triangle ID are set to -1,
  anchor vector to (0,0,0) -> the resulting data mark an invalid
  Interaction
*/
BTriMeshInteractionCpData::BTriMeshInteractionCpData() : 
  m_ap(),m_tid(-1),m_pid(-1)
{}

/*!
  Construct CpData from bonded triangle interaction - takes particle and triangle ID and particle anchor point

  \param BTI the bonded triangle interaction (reference)
*/
BTriMeshInteractionCpData::BTriMeshInteractionCpData(const BTriangleInteraction &BTI)
{
  m_pid=BTI.getPid();
  m_tid=BTI.getTid();
  m_ap=BTI.getAP();
}

/*!
  Set the data of an existing CpData object to those of a given 
  bonded triangle interaction.

  \param BTI the bonded triangle interaction (reference)
*/
void BTriMeshInteractionCpData::set(const BTriangleInteraction& BTI)
{
  m_pid=BTI.getPid();
  m_tid=BTI.getTid();
  m_ap=BTI.getAP();
}

/*!
  Write restartable CpData to output stream. The format is
  tid pid ap_x ap_y ap_z
  where ap_? are the vector components of the anchor point

  \warning The format isn't set in stone - don't rely on this for reading functions at the moment
*/ 
void BTriMeshInteractionCpData::saveCheckPointData(std::ostream& ost)
{
  const char delim = ' ';
  ost << m_pid << delim
      << m_tid << delim
      << m_ap;
}

/*
   Write visualization CpData to output stream.
*/
void BTriMeshInteractionCpData::saveSnapShotData(std::ostream& ost){
  const char delim = ' ';
  ost << m_pid << delim
      << m_tid << delim
      << m_ap;
}

/*!
  Load CpData for bonded triangle interaction from input stream

*/
void BTriMeshInteractionCpData::loadCheckPointData(std::istream& ist)
{}
