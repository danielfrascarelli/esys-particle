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

#include "RotElasticInteraction.h"
#include "console.h"

// --- TML includes ---
#include "tml/message/packed_message_interface.h"

#include <stdexcept>

CRotElasticIGP::CRotElasticIGP() : AIGParam(), m_kr(0.0), m_scaling(true)
{}

CRotElasticIGP::CRotElasticIGP(
  const  std::string& name,
  double kr,
  bool   scaling
) :
  AIGParam(name), m_kr(kr), m_scaling(scaling)
{}
  
CRotElasticInteraction::CRotElasticInteraction() : ARotPairInteraction()
{
  m_kr      = 0.0;
  m_nForce  = 0.0;
  m_force   = Vec3(0.0,0.0,0.0);
  m_scaling = true;
}

CRotElasticInteraction::CRotElasticInteraction(
  CRotParticle* p1,
  CRotParticle* p2,
  const CRotElasticIGP& param
) :
  ARotPairInteraction(p1,p2)
{
  double f=1.0; 
  m_scaling = param.m_scaling;
  // scale elastic param
  if (m_scaling) {
    if(!CParticle::getDo2dCalculations()) {
      f=0.5*(p1->getRad()+p2->getRad());
    }
  }
  m_kr      = f * param.m_kr;
  m_nForce  = 0.0;
  m_force   = Vec3(0.0,0.0,0.0);
  m_D       = p1->getPos()-p2->getPos();
}

Vec3 CRotElasticInteraction::getForce() const
{
  return m_force;
}

/*!
  Calculate free elastic forces.
*/
void CRotElasticInteraction::calcForces()
{ 
  const Vec3 rb = m_p1->getPos() - m_p2->getPos();
  const double rbNorm = rb.norm();
  const double r_0Norm = m_p1->getRad()+m_p2->getRad();  
  const Vec3 delta_r = (rbNorm - r_0Norm)*rb/rbNorm;

  if (rbNorm < r_0Norm) { // contact -> calculate forces
    m_force = m_kr*delta_r;
    m_nForce = -m_force.norm();

    Vec3 pos=m_p2->getPos()+(m_p2->getRad()/(m_p1->getRad()+m_p2->getRad()))*rb;
    m_p1->applyForce(-1.0*m_force,pos);
    m_p2->applyForce(m_force,pos);
    m_cpos=pos;
  } else { // no contact -> all forces are 0
    m_force = Vec3(0.0,0.0,0.0);
    m_nForce = 0.0;
  }
}

/*!
  Get the potential energy stored in the interaction.
*/
double CRotElasticInteraction::getPotentialEnergy() const
{
  return (m_kr!=0.0) ? 0.5*(m_nForce*m_nForce)/m_kr : 0.0;
}

/*!
  Get the particle member function which returns a scalar field of a given name.
 
  \param name the name of the field
*/
CRotElasticInteraction::ScalarFieldFunction CRotElasticInteraction::getScalarFieldFunction(const string& name)
{
  CRotElasticInteraction::ScalarFieldFunction sf;
  
  if (name=="potential_energy") {
    sf=&CRotElasticInteraction::getPotentialEnergy;
  } else if (name=="count") {
    sf=&CRotElasticInteraction::Count;
  } else {
    sf=NULL;
    cerr << "ERROR - invalid name for interaction scalar access function" << endl;
  }
  
  return sf;
}

/*!
  Get the particle member function which returns a checked scalar field of a given name.

  \param name the name of the field
*/
CRotElasticInteraction::CheckedScalarFieldFunction CRotElasticInteraction::getCheckedScalarFieldFunction(const string& name)
{
  CRotElasticInteraction::CheckedScalarFieldFunction sf;

  sf=NULL;
  cerr << "ERROR - invalid name for interaction vector access function" << endl;

  return sf;
}


/*!
  Get the particle member function which returns a vector field of a given name.
 
  \param name the name of the field
*/
CRotElasticInteraction::VectorFieldFunction CRotElasticInteraction::getVectorFieldFunction(const string& name)
{
  CRotElasticInteraction::VectorFieldFunction vf;
  
  if (name=="force"){
    vf=&CRotElasticInteraction::getForce;
  } else {
    vf=NULL;
    cerr << "ERROR - invalid name for interaction vector access function" << endl;
  }
  
  return vf;
}

/*!
  save restart data to ostream

  \param oStream the output stream
*/
void CRotElasticInteraction::saveRestartData(std::ostream &oStream)
{
  oStream << m_id[0] << " ";
  oStream << m_id[1] << " ";
  oStream << m_init  << " ";
  oStream << m_kr << " ";
  oStream << m_scaling << " ";
  oStream << m_D;
}

/*!
  load restart data from stream

  \param iStream the input stream
*/
void CRotElasticInteraction::loadRestartData(std::istream &iStream)
{
  iStream >> m_id[0];
  iStream >> m_id[1];
  iStream >> m_init ;
  iStream >> m_kr;
  iStream >> m_scaling;
  iStream >> m_D;
}

ostream& operator<<(ostream& ost,const CRotElasticInteraction& BI)
{
  ost << "[" << BI.m_p1->getID() << " - " << BI.m_p2->getID() << "]";
  return ost;
}
