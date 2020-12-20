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


#include "Model/RotPairInteraction.h"

ARotPairInteraction::~ARotPairInteraction()
{
}

ARotPairInteraction::ARotPairInteraction():AInteraction()
{
  m_p1=NULL;
  m_p2=NULL;
  m_id.clear();
  m_id.push_back(-1);
  m_id.push_back(-1);
}

ARotPairInteraction::ARotPairInteraction(CRotParticle* p1,CRotParticle* p2) : AInteraction()
{
  setPP(p1,p2);
  m_init=true;
}

void ARotPairInteraction::setPP(CRotParticle *p1,CRotParticle *p2)
{
  m_p1=p1;
  m_p2=p2;
  m_id.clear();
  m_id.push_back(m_p1->getID());
  m_id.push_back(m_p2->getID());
}

void ARotPairInteraction::setPP(const vector<CRotParticle*> pp)
{
  setPP(pp[0], pp[1]);
}

void ARotPairInteraction::checkIDs()
{
  if ((m_id[0]!=m_p1->getID())||(m_id[1]!=m_p2->getID())){
    std::cout
      << "inconsistent IDs : " << m_id[0] << "-" << m_id[1]
      << " vs. "
      << m_p1->getID() << "-" << m_p2->getID() << std::endl << std::flush;
  }
}

/*!
  check if any of the particles in the interaction fits tag & mask

  \param tag the tag
  \param mask the mask
*/
bool ARotPairInteraction::hasTag(int tag ,int mask) const
{
  int tag1=m_p1->getTag();
  int tag2=m_p2->getTag();
  return (((tag1 & mask)==(tag & mask))||((tag2 & mask)==(tag & mask)));
}
