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

#include "LinearDashpotInteraction.h"

// --- system includes ---
#include <iostream>

// --- member functions for interaction group parameter class ---

//! default constructor
CLinearDashpotIGP::CLinearDashpotIGP(): AIGParam(), m_damp(0.0), m_cutoff(1.0)
{}

/*!
  constructor with parameters

  \param name the name of the interaction group
  \param damp the damping coefficient
  \param cutoff the interaction range, relative to the sum of the particle radii
*/
CLinearDashpotIGP::CLinearDashpotIGP(const std::string& name,double damp, double cutoff)
  : AIGParam(name), m_damp(damp), m_cutoff(cutoff)
{}

// --- function of interaction class ---
CLinearDashpotInteraction::CLinearDashpotInteraction(CParticle* p1,CParticle* p2, const CLinearDashpotIGP& param) : APairInteraction(p1,p2)
{
  // calc. cross section - 2D / 3D
  double r_avg=0.5*(m_p1->getRad()+m_p2->getRad());
  if(CParticle::getDo2dCalculations()){
    m_cross_section=2.0*r_avg;
  } else {
    m_cross_section=M_PI*r_avg*r_avg;
  }
  // set parameters
  m_cutoff=param.m_cutoff;
  m_damp=param.m_damp;
  m_force=Vec3(0.0,0.0,0.0);
}


/*!
  calculate forces
*/
void CLinearDashpotInteraction::calcForces()
{
  // calc distance -> needed to check if particles interact
  Vec3 D=m_p1->getPos()-m_p2->getPos();
  double dist_sq=D*D;
  // calc cutoff distance -> could be cached
  double cut_dist=m_cutoff*(m_p1->getRad()+m_p2->getRad());
  if(dist_sq<(cut_dist*cut_dist)){
    // velocity difference
    Vec3 dvel=m_p1->getVel()-m_p2->getVel();
    // strain rate
    Vec3 eps_dot=dvel/sqrt(dist_sq);
    m_force=eps_dot*m_damp*m_cross_section;
    // apply force at particle centers
    m_p2->applyForce(m_force,m_p2->getPos());
    m_p1->applyForce(-1.0*m_force,m_p1->getPos()); 
  }
  m_cpos=(m_p1->getPos()+m_p2->getPos())*0.5;
}

/*!
  "field function" returning force currently exerted by interaction
*/
Vec3 CLinearDashpotInteraction::getForce() const
{
  return m_force;
}


/*!
  Get the particle member function which returns a scalar field of a given name.
 
  \param name the name of the field
*/
CLinearDashpotInteraction::ScalarFieldFunction CLinearDashpotInteraction::getScalarFieldFunction(const string& name)
{
  CLinearDashpotInteraction::ScalarFieldFunction sf;

  if (name=="count"){
    sf=&CLinearDashpotInteraction::Count;
  } else {
    sf=NULL;
    std::cerr << "ERROR - invalid name for interaction scalar  access function " << name << " in LinearDashpotInteraction" << std::endl;
  }

  return sf;
}

/*!
  Get the particle member function which returns a vector field of a given name.
 
  \param name the name of the field
*/
CLinearDashpotInteraction::VectorFieldFunction CLinearDashpotInteraction::getVectorFieldFunction(const string& name)
{
  CLinearDashpotInteraction::VectorFieldFunction vf;

  if (name=="force"){
    vf=&CLinearDashpotInteraction::getForce;
  } else {
    vf=NULL;
    std::cerr << "ERROR - invalid name for interaction vector access function " << name << " in LinearDashpotInteraction"  << endl;
  }
  
  return vf;
}
 
/*!
  dummy
*/
CLinearDashpotInteraction::CheckedScalarFieldFunction CLinearDashpotInteraction::getCheckedScalarFieldFunction(const string& name)
{
  CLinearDashpotInteraction::CheckedScalarFieldFunction csf;

  csf=NULL;
  cerr << "ERROR - invalid name for interaction vector access function " << name << " in LinearDashpotInteraction"  << endl;
  
  return csf;
}
