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


#include "Parallel/mpibuf.h"
#include "Model/Interaction.h"

#include <iostream>

AInteraction::AInteraction(): m_id(), m_iid(-1), m_init(false)
{
}

AInteraction::~AInteraction()
{
}

vector<int> AInteraction::getAllID() const
{
  return m_id;
}

bool AInteraction::initialized() const
{
  return m_init;
}

APairInteraction::~APairInteraction()
{
}

APairInteraction::APairInteraction():AInteraction()
{
  m_p1=NULL;
  m_p2=NULL;
  m_id.clear();
  m_id.push_back(-1);
  m_id.push_back(-1);
}

APairInteraction::APairInteraction(CParticle* p1, CParticle* p2)
  : AInteraction()
{
  // force the particle with the lower ID to be p1
  if(p1->getID()<p2->getID()){
    m_p1=p1;
    m_p2=p2;
  } else {
    m_p1=p2;
    m_p2=p1;
  }
  m_id.clear();
  m_id.push_back(m_p1->getID());
  m_id.push_back(m_p2->getID());
  m_init=true;
}

void APairInteraction::setPP(CParticle *p1,CParticle *p2)
{
  m_p1=p1;
  m_p2=p2;
  m_id.clear();
  m_id.push_back(m_p1->getID());
  m_id.push_back(m_p2->getID());
}

void APairInteraction::checkIDs()
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
bool APairInteraction::hasTag(int tag ,int mask) const
{
  int tag1=m_p1->getTag();
  int tag2=m_p2->getTag();
  return (((tag1 & mask)==(tag & mask))||((tag2 & mask)==(tag & mask)));
}
