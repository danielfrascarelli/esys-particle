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

#include "EEdgeInteraction.h"
#include "Foundation/vec3.h"

// -- STL includes --
#include <utility>
using std::pair;

/*!
  default constructor
*/
EEdgeInteraction::EEdgeInteraction()
{
  m_p=NULL;
  m_edge=NULL;
  m_k=0.0;
  m_inner_flag=false;
}

/*!
  constructor with parameters

  \param p
  \param e
  \param param
  \param iflag
*/
EEdgeInteraction::EEdgeInteraction(CParticle* p,Edge* e,ETriMeshIP param ,bool iflag)
{
  m_p=p;
  m_edge=e;
  // scale elastic param
  double f=1.0;
  if(!CParticle::getDo2dCalculations()){
    f*=3.141592654*this->m_p->getRad();
  }
  m_k=f*param.k;
  m_inner_flag=iflag;
}

/*!
  destructor
*/
EEdgeInteraction::~EEdgeInteraction()
{}


/*!
  calculate & apply forces
*/
void EEdgeInteraction::calcForces()
{
  Vec3 ppos=m_p->getPos();
  if(m_edge->isValidContact(ppos)){ // if no contact to adjacent triangles
    pair<bool,double> dist=m_edge->dist(ppos);
    if(dist.first && (dist.second<m_p->getRad())){
      Vec3 force=m_k*(m_p->getRad()-dist.second)*m_edge->getDirectionFromPoint(m_p->getPos());
      Vec3 pos=m_p->getPos()-dist.second*m_edge->getDirectionFromPoint(m_p->getPos());
      m_p->applyForce(force,pos);
      if(m_inner_flag) m_edge->applyForce(-1.0*force);
    } 
  } 
}
