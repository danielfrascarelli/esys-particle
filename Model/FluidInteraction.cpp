/////////////////////////////////////////////////////////////
//                                                         //
// Copyright (c) 2003-2014 by The University of Queensland //
// Centre for Geoscience Computing                         //
// http://earth.uq.edu.au/centre-geoscience-computing      //
//                                                         //
// Primary Business: Brisbane, Queensland, Australia       //
// Licensed under the Open Software License version 3.0    //
// http://www.apache.org/licenses/LICENSE-2.0          //
//                                                         //
/////////////////////////////////////////////////////////////


#include "Model/FluidInteraction.h"


CFluidInteraction::CFluidInteraction (CParticle* p,CFluidCell* cell)
{
  m_p=p;
  m_cell=cell;
}


/*!
Calculate the fluid force.
*/
void CFluidInteraction::calcForces()
{
  double Mu=m_cell->getMu();
  double Phi=m_cell->geteffPhi();
  double D=m_cell->getD();
  Vec3 vrel=m_cell->getVf()-m_p->getVel();
  //drag force on unit volume
  Vec3 drag;
  if(D==0){drag=Vec3(0,0,0);}
  else{drag=(150.0*Mu*(1.0-Phi)/Phi/D/D+1.75*1000.0*vrel.norm()/D)*vrel;};
  //buoyancy force on unit volume
  Vec3 buoyancy=m_cell->getPg();
  //force applied on particle
  double r=m_p->getRad();
  double Volume=4.0/3.0*3.1415926*r*r*r;
  m_drag=drag*Volume;
  m_buoyancy=buoyancy*Volume;
  m_force=m_drag+m_buoyancy;
  m_p->applyForce(m_force,m_p->getPos());
}


/*!
  Get the fluidinteraction member function which returns a scalar field of a given name.

  \param name the name of the field
*/
CFluidInteraction::ScalarFieldFunction CFluidInteraction::getScalarFieldFunction(const string& name)
{
  CFluidInteraction::ScalarFieldFunction sf;
  if(name=="ID"){
    sf=&CFluidInteraction::getParticleID;
  } else if(name=="Drag_abs"){
    sf=&CFluidInteraction::getVbsDrag;
  } else if(name=="Buoyancy_abs"){
    sf=&CFluidInteraction::getVbsBuoyancy;
  } else if(name=="Force_abs"){
    sf=&CFluidInteraction::getVbsForce;
  } else {
    sf=NULL;
    cerr << "ERROR - invalid name for fluid interaction scalar  access function" << endl;
  }

  return sf;
}


/*!
  Get the fluid interaction member function which returns a vector field of a given name.

  \param name the name of the field
*/
CFluidInteraction::VectorFieldFunction CFluidInteraction::getVectorFieldFunction(const string& name)
{
  CFluidInteraction::VectorFieldFunction sf;

  if(name=="Position"){
    sf=&CFluidInteraction::getParticlePos;
  } else if(name=="Drag"){
    sf=&CFluidInteraction::getDrag;
  } else if(name=="Buoyancy"){
    sf=&CFluidInteraction::getBuoyancy;
  } else if(name=="Force"){
    sf=&CFluidInteraction::getForce;
  } else {
    sf=NULL;
    cerr << "ERROR - invalid name for fluid interaction vector access function" << endl;
  }

  return sf;
}


