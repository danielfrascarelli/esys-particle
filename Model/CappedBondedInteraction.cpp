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

#include "Model/CappedBondedInteraction.h"
#include "Foundation/console.h"
#include "tml/message/packed_message_interface.h"

#include <stdexcept>

CCappedBondedIGP::CCappedBondedIGP() : CBondedIGP(),m_force_limit(0.0)
{
}
  
CCappedBondedIGP::CCappedBondedIGP(const std::string &name, int tag, double normalK, double breakDistance,double maxforce)
  : CBondedIGP(name,tag,normalK,breakDistance),
    m_force_limit(maxforce)
{
}

CCappedBondedInteraction::CCappedBondedInteraction() : CBondedInteraction()
{
  m_force_limit=0.0;
}

/*!
  just do the APairInteraction part of the constructor - not to be used directly, 
  only by derived class -> therefore protected
*/
CCappedBondedInteraction::CCappedBondedInteraction(CParticle* p1,CParticle* p2)
  : CBondedInteraction(p1, p2)
{}

CCappedBondedInteraction::CCappedBondedInteraction(
  CParticle* p1,
  CParticle* p2,
  const CCappedBondedIGP& param
)
  : CBondedInteraction(p1, p2)
{
    m_k=param.k;
  m_break=(m_p1->getRad()+m_p2->getRad())*param.rbreak;
  m_r0=p1->getRad()+p2->getRad();
  // setup initial distance in case we want to break before first calcForces()
  Vec3 D=p1->getPos()-p2->getPos();
  m_dist=sqrt(D*D);
  m_force=Vec3(0.0,0.0,0.0);
  m_tag = param.tag;
  m_force_limit=param.m_force_limit;
}

CCappedBondedInteraction::~CCappedBondedInteraction()
{
  //
  // if (m_p1!=NULL) m_p1->setFlag();
  // if (m_p2!=NULL) m_p2->setFlag();
}


/*!
  Calculate bonded elastic forces. 21 Flops
*/
void CCappedBondedInteraction::calcForces()
{
  //  console.XDebug() 
  //    << "bonded interaction: [" << m_p1->getID() << " - "
  //    << m_p2->getID() << "]" << m_p1->getPos()
  //    << m_p2->getPos() << "\n"; 
  const Vec3 D=m_p1->getPos()-m_p2->getPos();
  // get distance
  m_dist=sqrt(D*D);

  // calc relative deformation
  double rdef=(m_dist-m_r0)/m_dist;

  if(fabs(rdef)>m_force_limit){
    rdef=(rdef/fabs(rdef))*m_force_limit;
  }
  m_force = D*(m_k*rdef);

  const Vec3 pos=m_p2->getPos()+(m_p2->getRad()/m_r0)*D;

  m_p2->applyForce(m_force,pos);
  m_p1->applyForce(-1.0*m_force,pos);
  m_cpos=pos;
}

/*!
  Get the particle member function which returns a scalar field of a given name.

  \param name the name of the field 
*/
CCappedBondedInteraction::ScalarFieldFunction CCappedBondedInteraction::getScalarFieldFunction(const string& name)
{
  CCappedBondedInteraction::ScalarFieldFunction sf;

  if (name=="potential_energy"){
    sf=&CBondedInteraction::getPotentialEnergy;
  } else if (name=="count"){
    sf=&CBondedInteraction::Count;
  } else if (name=="strain"){
    sf=&CBondedInteraction::getStrain;;
  } else {
    sf=NULL;
    cerr << "ERROR - invalid name for interaction scalar  access function" << endl; 
  }
  
  return sf;
}

/*!
  Get the particle member function which returns a vector field of a given name.

  \param name the name of the field 
*/
CCappedBondedInteraction::VectorFieldFunction CCappedBondedInteraction::getVectorFieldFunction(const string& name)
{
  CCappedBondedInteraction::VectorFieldFunction sf;

  if (name=="force"){
    sf=&CBondedInteraction::getForce;
  } else {
    sf=NULL;
    cerr << "ERROR - invalid name for interaction vector  access function" << endl; 
  }

  return sf;
}

/*!
  Get the particle member function which returns a checked scalar field of a given name.

  \param name the name of the field
*/
CCappedBondedInteraction::CheckedScalarFieldFunction CCappedBondedInteraction::getCheckedScalarFieldFunction(const string& name)
{
	CCappedBondedInteraction::CheckedScalarFieldFunction sf;

	sf=NULL;
	cerr << "ERROR - invalid name for interaction scalar  access function" << endl;

	return sf;
}

/*!
  Pack a CBondedInteraction into a TML packed message

  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::pack<CCappedBondedInteraction>(const CCappedBondedInteraction& I)
{
  append(I.m_k);
  append(I.m_r0);
  append(I.m_dist);
  append(I.m_break);
  append(I.m_id[0]);
  append(I.m_id[1]);
  append(I.getTag());
  append(I.m_force_limit);
}

/*!
  Unpack a CBondedInteraction from a TML packed message

  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::unpack<CCappedBondedInteraction>(CCappedBondedInteraction& I)
{
  I.m_k=pop_double();
  I.m_r0=pop_double();
  I.m_dist=pop_double();
  I.m_break=pop_double();
  I.m_id.erase(I.m_id.begin(),I.m_id.end());
  I.m_id.push_back(pop_int());
  I.m_id.push_back(pop_int());
  I.setTag(pop_int());
  I.m_force_limit=pop_double();
}

