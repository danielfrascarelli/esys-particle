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

#include "HertzianElasticInteraction.h"

// --- system includes ---
#include <iostream>

// --- member functions for interaction group parameter class ---

//! default constructor
CHertzianElasticIGP::CHertzianElasticIGP(): AIGParam(), m_E(0.0), m_nu(0.0)
{}

/*!
  constructor with parameters

  \param name the name of the interaction group
  \param k the stiffness constant
*/
CHertzianElasticIGP::CHertzianElasticIGP(const std::string& name,double E,double nu)
  : AIGParam(name), m_E(E), m_nu(nu)
{}

// --- function of interaction class ---
CHertzianElasticInteraction::CHertzianElasticInteraction(CParticle* p1,CParticle* p2, const CHertzianElasticIGP& param) : APairInteraction(p1,p2)
{
  m_E=param.m_E;
  m_nu=param.m_nu;
  m_force=Vec3(0.0,0.0,0.0);
  m_dn=0.0;
}


/*!
  calculate forces
*/
void CHertzianElasticInteraction::calcForces()
{
  Vec3 D=m_p1->getPos()-m_p2->getPos();
  double dist=D*D;
  double eq_dist=m_p1->getRad()+m_p2->getRad();
  //std::cerr << "d: " << dist << " / " << eq_dist*eq_dist << std::endl;
  if(dist<(eq_dist*eq_dist)){
    double R_ij=1.0/(1.0/m_p1->getRad()+1.0/m_p2->getRad());
    dist=sqrt(dist);
    m_dn=eq_dist-dist;
    Vec3 dir=D.unit();
    m_force=dir*(m_E*sqrt(R_ij))/(2.0*(1.0-m_nu*m_nu))*pow(m_dn,1.5);
    //std::cerr << "HF: " <<  dir << " " << R_ij << " " << m_dn << " " << m_force << std::endl; 
    Vec3 pos=m_p2->getPos()+(m_p2->getRad()/eq_dist)*D;
    m_p1->applyForce(m_force,pos);
    m_p2->applyForce(-1.0*m_force,pos); 
  } else {
    m_force=Vec3(0.0,0.0,0.0);
    m_dn=0.0;
  }
}

/*!
  "field function" returning force currently exerted by interaction
*/
Vec3 CHertzianElasticInteraction::getForce() const
{
  return m_force;
}

/*!
  "field function" returning potential energy currently stored in interaction
*/
double CHertzianElasticInteraction::getPotentialEnergy() const
{
  return 0.4*m_force.norm()*m_dn;
}


/*!
  Get the particle member function which returns a scalar field of a given name.
 
  \param name the name of the field
*/
CHertzianElasticInteraction::ScalarFieldFunction CHertzianElasticInteraction::getScalarFieldFunction(const string& name)
{
  CHertzianElasticInteraction::ScalarFieldFunction sf;

  if (name=="potential_energy"){
    sf=&CHertzianElasticInteraction::getPotentialEnergy;
  } else if (name=="count"){
    sf=&CHertzianElasticInteraction::Count;
  } else {
    sf=NULL;
    std::cerr << "ERROR - invalid name for interaction scalar  access function" << std::endl;
  }
  
  return sf;
}

/*!
  Get the particle member function which returns a vector field of a given name.
 
  \param name the name of the field
*/
CHertzianElasticInteraction::VectorFieldFunction CHertzianElasticInteraction::getVectorFieldFunction(const string& name)
{
  CHertzianElasticInteraction::VectorFieldFunction vf;

  if (name=="force"){
    vf=&CHertzianElasticInteraction::getForce;
  } else {
    vf=NULL;
    cerr << "ERROR - invalid name for interaction vector access function" << endl;
  }
  
  return vf;
}
 
/*!
  dummy
*/
CHertzianElasticInteraction::CheckedScalarFieldFunction CHertzianElasticInteraction::getCheckedScalarFieldFunction(const string& name)
{
  CHertzianElasticInteraction::CheckedScalarFieldFunction csf;

  csf=NULL;
  cerr << "ERROR - invalid name for interaction vector access function" << endl;
  
  return csf;
}
