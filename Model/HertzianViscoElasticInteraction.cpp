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

/*
HertzianViscoElasticInteraction.cpp:
  Written by Laura Heredia and Pablo Richeri, 2009.
*/

#include "HertzianViscoElasticInteraction.h"

// --- system includes ---
#include <iostream>

// --- member functions for interaction group parameter class ---

//! default constructor
CHertzianViscoElasticIGP::CHertzianViscoElasticIGP()
  : AIGParam(), m_A(0.0), m_E(0.0), m_nu(0.0)
{
}

/*!
  constructor with parameters

  \param name the name of the interaction group
  \param k the stiffness constant
*/
CHertzianViscoElasticIGP::CHertzianViscoElasticIGP(
  const std::string& name,
  double A,
  double E,
  double nu
)
  : AIGParam(name), m_A(A), m_E(E), m_nu(nu)
{
}

// --- function of interaction class ---
CHertzianViscoElasticInteraction::CHertzianViscoElasticInteraction(
  CParticle* p1,
  CParticle* p2,
  const CHertzianViscoElasticIGP& param
)
  : APairInteraction(p1,p2)
{
  m_A=param.m_A;
  m_E=param.m_E;
  m_nu=param.m_nu;
  m_force=Vec3(0.0,0.0,0.0);
  m_dn=0.0;
}


/*!
  calculate forces
*/
void CHertzianViscoElasticInteraction::calcForces()
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

    //Calculate d m_dn / dt
    double ex=D.X()/dist;
    double ey=D.Y()/dist;
    double ez=D.Z()/dist;
    double dvx=m_p1->getVel().X()-m_p2->getVel().X();
    double dvy=m_p1->getVel().Y()-m_p2->getVel().Y();
    double dvz=m_p1->getVel().Z()-m_p2->getVel().Z();
    double m_dn_dot=-(ex*dvx+ey*dvy+ez*dvz);
    
    double norm_m_force =
      (2.0*m_E*sqrt(R_ij)) /
      (3.0*(1.0-m_nu*m_nu)) *
      (pow(m_dn,1.5)+m_A*sqrt(m_dn)*m_dn_dot);
    m_force=norm_m_force<0?Vec3(0.0,0.0,0.0):dir*norm_m_force;
    //std::cerr << "HF: " <<  dir << " " << R_ij << " " << m_dn << " "
    //  << m_force << std::endl; 
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
Vec3 CHertzianViscoElasticInteraction::getForce() const
{
  return m_force;
}

/*!
  "field function" returning potential energy currently stored in interaction
*/
double CHertzianViscoElasticInteraction::getPotentialEnergy() const
{
  return 0.4*m_force.norm()*m_dn;
}


/*!
  Get the particle member function which returns a scalar field of a given name.
 
  \param name the name of the field
*/
CHertzianViscoElasticInteraction::ScalarFieldFunction
CHertzianViscoElasticInteraction::getScalarFieldFunction(const string& name)
{
  CHertzianViscoElasticInteraction::ScalarFieldFunction sf;

  if (name=="potential_energy"){
    sf=&CHertzianViscoElasticInteraction::getPotentialEnergy;
  } else if (name=="count"){
    sf=&CHertzianViscoElasticInteraction::Count;
  } else {
    sf=NULL;
    std::cerr
      << "ERROR - invalid name for interaction scalar  access function"
      << std::endl;
  }
  
  return sf;
}

/*!
  Get the particle member function which returns a vector field of a given name.
 
  \param name the name of the field
*/
CHertzianViscoElasticInteraction::VectorFieldFunction
CHertzianViscoElasticInteraction::getVectorFieldFunction(const string& name)
{
  CHertzianViscoElasticInteraction::VectorFieldFunction vf;

  if (name=="force"){
    vf=&CHertzianViscoElasticInteraction::getForce;
  } else {
    vf=NULL;
    cerr
      << "ERROR - invalid name for interaction vector access function"
      << endl;
  }
  
  return vf;
}
 
/*!
  dummy
*/
CHertzianViscoElasticInteraction::CheckedScalarFieldFunction
CHertzianViscoElasticInteraction::getCheckedScalarFieldFunction(
  const string& name
)
{
  CHertzianViscoElasticInteraction::CheckedScalarFieldFunction csf;

  csf=NULL;
  cerr
    << "ERROR - invalid name for interaction vector access function"
    << endl;
  
  return csf;
}
