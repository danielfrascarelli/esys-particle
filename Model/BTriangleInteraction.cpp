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

#include "BTriangleInteraction.h"
#include "tml/message/packed_message_interface.h"

/*!
  default constructor
*/
BTriangleInteraction::BTriangleInteraction()
{
  m_p=NULL;
  m_t=NULL;
  m_tid=-1;
  m_pid=-1;
}


/*!
  constructor with parameters

  \param p a pointer to the particle
  \param t a pointer to the triangle
  \param param the interaction parameters
  \param iflag 
*/
BTriangleInteraction::BTriangleInteraction(CParticle* p,Triangle* t,BTriMeshIP param,bool iflag)
{
  m_p=p;
  m_t=t;
  m_k=param.k;
  m_break=param.brk*m_p->getRad();
  m_inner_flag=iflag;
  // setup anchor point
  m_ap=m_t->toLocal(m_p->getPos());
  m_dist=0.0; // inital distance is always 0.0 !
  m_tid=m_t->getID();
  m_pid=m_p->getID();
}


/*!
  destructor
*/
BTriangleInteraction::~BTriangleInteraction()
{}

/*!
  calculate & apply forces
*/
void BTriangleInteraction::calcForces()
{
  // transform anchor point to world coords
  Vec3 ap_global=m_t->toGlobal(m_ap);
  // get dist between anchor and particle
  const Vec3 D=ap_global-m_p->getPos();
  m_dist=sqrt(D*D);
  // calc force
  Vec3 force=D*m_k;
  Vec3 pos=m_p->getPos();
  // apply force
  m_p->applyForce(force,pos);
  if(m_inner_flag) m_t->applyForce(-1.0*force);
}

/*!
  return if the interaction is broken, i.e. the distance between 
  particle and anchor point exceeds breaking distance, i.e. relative
  breaking distance x particle readius
*/
bool BTriangleInteraction::broken()
{
  return (m_dist>m_break);
}


/*!
  returns the projection of the anchor point on the triangle
*/
Vec3 BTriangleInteraction::getAP() const
{ 
  Vec3 app=Vec3(m_ap.X(),m_ap.Y(),0.0); // move point onto triangle
  Vec3 ap_global=m_t->toGlobal(app); // transfrom point to global coordinates

  return ap_global;
}

/*!
  Pack a BTriangleInteraction into a TML packed message

  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::pack<BTriangleInteraction>(const BTriangleInteraction& I)
{
  append(I.m_k);
  append(I.m_dist);
  append(I.m_break);
  append(I.m_ap);
  append(I.getTid());
  append(I.getPid());
}

/*!
  Unpack a BTriangleInteraction from a TML packed message

  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::unpack<BTriangleInteraction>(BTriangleInteraction& I)
{
  I.m_k=pop_double();
  I.m_dist=pop_double();
  I.m_break=pop_double();
  I.m_ap=pop_vec3();
  I.m_tid=pop_int();
  I.m_pid=pop_int();
}

