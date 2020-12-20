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

#include "Model/BondedInteraction.h"
#include "Model/BondedInteractionCpData.h"
#include "Foundation/console.h"
#include "tml/message/packed_message_interface.h"

#include <stdexcept>

CBondedIGP::CBondedIGP() : AIGParam(), k(0.0), rbreak(0.0), tag(0)
{
}
  
CBondedIGP::CBondedIGP(const std::string &name, int tag, double normalK, double breakDistance, bool scaling)
  : AIGParam(name),
    k(normalK),
    rbreak(breakDistance),
    tag(tag),
    m_scaling(scaling)
{
}

CBondedInteraction::CBondedInteraction() : APairInteraction()
{
  m_k=0;
  m_break=0;
  m_r0=0;
  m_dist=0;
  m_force=Vec3(0.0,0.0,0.0);
  m_tag = -1;
  m_scaling=true;
}

/*!
  just do the APairInteraction part of the constructor - not to be used directly, 
  only by derived class -> therefore protected
*/
CBondedInteraction::CBondedInteraction(CParticle* p1,CParticle* p2)
  : APairInteraction(p1, p2)
{}

CBondedInteraction::CBondedInteraction(
  CParticle* p1,
  CParticle* p2,
  const CBondedIGP& param
)
  : APairInteraction(p1, p2)
{
  double effR=1.0;
  double effA=1.0; 
  double effL=1.0; // effective radius, cross section  and length of the bond for scaling 

  // equilibrium distance
  m_r0=p1->getRad()+p2->getRad();
  
  // scale elastic param
  m_scaling = param.m_scaling;
  if (m_scaling) {
    if(!CParticle::getDo2dCalculations()){
      effR=0.5*m_r0;
    }
    effL=m_r0;
    effA = effR * effR;
  }
  m_k=param.k*effA/effL;  
  m_break=m_r0*param.rbreak;
  // setup initial distance in case we want to break before first calcForces()
  Vec3 D=p1->getPos()-p2->getPos();
  m_dist=sqrt(D*D);
  m_force=Vec3(0.0,0.0,0.0);
  m_tag = param.tag;
  m_scaling = param.m_scaling;
}

CBondedInteraction::~CBondedInteraction()
{
}

bool CBondedInteraction::broken()
{
  bool res;

  if((m_dist-m_r0)>m_break){
//     console.Debug() << "bond broken" << "\n" 
// 		    << "ids : " << m_p1->getID() <<  " " << m_p2->getID() << "\n"
// 		    << "positions : " << m_p1->getPos() <<  m_p2->getPos() << "\n"
// 		    << "dist : " << m_dist << "\n" 
// 		    << "break : " << m_break << "\n" ;
    res=true;
    if (m_p1!=NULL) m_p1->setFlag();
    if (m_p2!=NULL) m_p2->setFlag();
  } else {
    res=false;
  }
  return(res);
}

double CBondedInteraction::getCriterion() const
{
  return m_dist - m_break;
}

/*!
  Set breaking distance to sum of the radii multiplied with given "relative breaking distance"

  \param rel_break the relative breaking distance
*/
void CBondedInteraction::setBreak(double rel_break)
{
  m_break=(m_p1->getRad()+m_p2->getRad())*rel_break;
}

/*!
  Calculate bonded elastic forces. 21 Flops
*/
void CBondedInteraction::calcForces()
{
  const Vec3 D=m_p1->getPos()-m_p2->getPos();
  m_dist=sqrt(D*D);

  m_force = D*(m_k*(m_dist-m_r0)/m_dist);
 //  console.XDebug() 
//     << "bonded interaction: [" << m_p1->getID() << " - "
//     << m_p2->getID() << "]" << m_p1->getPos()
//     << m_p2->getPos() << " " << m_force << "\n"; 
  
  const Vec3 pos=m_p2->getPos()+(m_p2->getRad()/m_r0)*D;

  m_p2->applyForce(m_force,pos);
  m_p1->applyForce(-1.0*m_force,pos);
  m_cpos=pos;
}

/*!
  get the potential energy stored in the interaction
*/
double CBondedInteraction::getPotentialEnergy() const
{
  double e_pot_norm=0.5*m_force*m_force/m_k;

  return e_pot_norm;
}


/*!
  get strain - compression positive
*/
double CBondedInteraction::getStrain() const
{
  double strain=(m_r0-m_dist)/m_r0;

  return strain;
}

/*!
  get force - points to p1 on extension, to p2 on compression
*/
Vec3 CBondedInteraction::getForce() const
{
  return m_force;
}

/*!
  Get the particle member function which returns a scalar field of a given name.

  \param name the name of the field 
*/
CBondedInteraction::ScalarFieldFunction CBondedInteraction::getScalarFieldFunction(const string& name)
{
  CBondedInteraction::ScalarFieldFunction sf;

  if (name=="potential_energy"){
    sf=&CBondedInteraction::getPotentialEnergy;
  } else if (name=="count"){
    sf=&CBondedInteraction::Count;
  } else if (name=="strain"){
    sf=&CBondedInteraction::getStrain;;
  } else if (name=="breaking_criterion"){
    sf=&CBondedInteraction::getCriterion;
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
CBondedInteraction::VectorFieldFunction CBondedInteraction::getVectorFieldFunction(const string& name)
{
  CBondedInteraction::VectorFieldFunction vf;

  if (name=="force"){
    vf=&CBondedInteraction::getForce;
  } else {
    vf=NULL;
    cerr << "ERROR - invalid name for interaction vector  access function" << endl; 
  }

  return vf;
}

/*!
  Get the particle member function which returns a checked scalar field of a given name.

  \param name the name of the field
*/
CBondedInteraction::CheckedScalarFieldFunction CBondedInteraction::getCheckedScalarFieldFunction(const string& name)
{
	CBondedInteraction::CheckedScalarFieldFunction sf;

	sf=NULL;
	cerr << "ERROR - invalid name for interaction scalar  access function" << endl;

	return sf;
}

/*!
  Pack a CBondedInteraction into a TML packed message

  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::pack<CBondedInteraction>(const CBondedInteraction& I)
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
void TML_PackedMessageInterface::unpack<CBondedInteraction>(CBondedInteraction& I)
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

/*!
  Save snapshot data (non-restartable, viz/postprocessing only) to an
  output stream.

  \param  oStream the output stream
*/
void CBondedInteraction::saveCheckPointData(std::ostream &oStream)
{
  BondedInteractionCpData(*this).saveCheckPointData(oStream);
}

/*!
  save restart data to ostream

  \param oStream the output stream
*/
void CBondedInteraction::saveRestartData(std::ostream &oStream)
{
  oStream << m_k << " ";
  oStream << m_r0 << " ";
  oStream << m_dist  << " ";
  oStream << m_break << " ";
  oStream << m_scaling << " ";
  oStream << m_id[0] << " ";
  oStream << m_id[1] << " ";
  oStream << getTag() << std::endl;
}

/*!
  load restart data from stream

  \param iStream the input stream
*/
void CBondedInteraction::loadRestartData(std::istream &iStream)
{
  int tag;
  iStream >> m_k ;
  iStream >> m_r0 ;
  iStream >> m_dist  ;
  iStream >> m_break ;
  iStream >> m_scaling ;
  iStream >> m_id[0] ;
  iStream >> m_id[1] ;
  iStream >> tag;
  setTag(tag);
}

ostream& operator<<(ostream& ost,const CBondedInteraction& BI)
{
  ost << "[" << BI.m_p1->getID() << " - " << BI.m_p2->getID() << "]";
  return ost;
}
