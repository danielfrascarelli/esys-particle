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

#include "Model/ShortBondedInteraction.h"
#include "tml/message/packed_message_interface.h"
#include "Foundation/console.h"

// --- STL includes ---
#include <stdexcept>

using std::runtime_error;

CShortBondedInteraction::CShortBondedInteraction() 
{}

/*!
  Construct valid short bonded interaction. The equilibrium distance is calculated from the 
  initial distance of the two particles. 

  \param p1 pointer to 1st particle
  \param p2 pointer to 2nd particles
  \param param the interaction parameters
*/
CShortBondedInteraction::CShortBondedInteraction(CParticle* p1,CParticle* p2,const CBondedIGP& param)
  : CBondedInteraction(p1,p2)
{
  m_k=param.k;
  m_r0=(p1->getPos()-p2->getPos()).norm();;
  m_break=m_r0*param.rbreak;
  m_dist=m_r0;;
  m_force=Vec3(0.0,0.0,0.0);
  m_tag = param.tag;
}

CShortBondedInteraction::~CShortBondedInteraction()
{}

CShortBondedInteraction::ScalarFieldFunction CShortBondedInteraction::getScalarFieldFunction(const string& name)
{
  CShortBondedInteraction::ScalarFieldFunction sf;

  if (name=="potential_energy"){
    sf=&CShortBondedInteraction::getPotentialEnergy;
  } else if (name=="count"){
    sf=&CShortBondedInteraction::Count;
  } else if (name=="strain"){
    sf=&CShortBondedInteraction::getStrain;;
  } else {
    sf=NULL;
    cerr << "ERROR - invalid name for interaction scalar  access function" << endl; 
  }
  
  return sf;

}

CShortBondedInteraction::CheckedScalarFieldFunction CShortBondedInteraction::getCheckedScalarFieldFunction(const string& name)
{
  CShortBondedInteraction::CheckedScalarFieldFunction sf;

  sf=NULL;
  cerr << "ERROR - invalid name for interaction scalar  access function" << endl;
  
  return sf;
}

CShortBondedInteraction::VectorFieldFunction CShortBondedInteraction::getVectorFieldFunction(const string& name)
{
  CShortBondedInteraction::VectorFieldFunction sf;

  if (name=="force"){
    sf=&CShortBondedInteraction::getForce;
  } else {
    sf=NULL;
    cerr << "ERROR - invalid name for interaction vector  access function" << endl; 
  }

  return sf;
}

void CShortBondedInteraction::saveCheckPointData(std::ostream &oStream)
{
  ShortBondedInteractionCpData(*this).saveCheckPointData(oStream);
}

void CShortBondedInteraction::loadCheckPointData(std::istream &iStream)
{
  throw runtime_error("CShortBondedInteraction::loadCheckPointData not implemented.");
}

ostream& operator<<(ostream& ost,const CShortBondedInteraction& BI)
{
  ost << CBondedInteraction(BI);
  return ost;
}

/*!
  Pack a CBondedInteraction into a TML packed message

  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::pack<CShortBondedInteraction>(const CShortBondedInteraction& I)
{
  append(I.m_k);
  append(I.m_r0);
  append(I.m_dist);
  append(I.m_break);
  append(I.m_id[0]);
  append(I.m_id[1]);
  append(I.getTag());
}

/*!
  Unpack a CBondedInteraction from a TML packed message

  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::unpack<CShortBondedInteraction>(CShortBondedInteraction& I)
{
  I.m_k=pop_double();
  I.m_r0=pop_double();
  I.m_dist=pop_double();
  I.m_break=pop_double();
  I.m_id.erase(I.m_id.begin(),I.m_id.end());
  I.m_id.push_back(pop_int());
  I.m_id.push_back(pop_int());
  I.setTag(pop_int());
}
