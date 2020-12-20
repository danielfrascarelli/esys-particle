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

#include "BCorner2DInteraction.h"

// -- Project includes --
#include "Foundation/vec3.h"
#include "tml/message/packed_message_interface.h"
#include "console.h"

/*!
  default constructor
*/
BCorner2DInteraction::BCorner2DInteraction()
{
  m_p=NULL;
  m_corner=NULL;
  m_pid=-1;
  m_cid=-1;
}

/*!
  constructor with parameters

  \param p a pointer to the particle
  \param c a pointer to the corner
  \param param the interaction parameters
  \param iflag 
*/
BCorner2DInteraction::BCorner2DInteraction(CParticle* p,Corner2D* c,BMesh2DIP param,bool iflag)
{
  m_p=p;
  m_corner=c;
  m_k=param.k;
  m_break=param.brk*m_p->getRad();
  // setup anchor point coefficients
  int m_ne=m_corner->getNEdges();
  if (m_ne==1){ // single edge case
    console.Critical() << "Signle Edge Case not implemented\n";
  } else if (m_ne==2){ // two edge (normal) case
    Vec3 n1=m_corner->getEdgeNormal(1);
    Vec3 n2=m_corner->getEdgeNormal(2);
    Vec3 p=m_p->getPos()-m_corner->getPos();
    k1=(n2.Y()*p.X()-n2.X()*p.Y())/(n1.X()*n2.Y()-n1.Y()*n2.X());
    k2=(n1.Y()*p.X()-n1.X()*p.Y())/(n2.X()*n1.Y()-n2.Y()*n1.X());

    // check
    Vec3 check=k1*n1+k2*n2;
    console.XDebug() << "BCorner2DInteraction check: " << check-p << "\n";
//     cout << "BCorner2DInteraction check: n1,n2,p,k1,k2 [" << n1 << "] [" << n2 << "] [" << p << "] , " << k1  << " , " << k2 << " [" << check-p<< "]\n";
  } else {
    console.Critical() << "ERROR: Corner appears to have 0 Edges\n";
  }
  m_dist=0.0; // inital distance is always 0.0 !
  m_pid=m_p->getID();
  m_cid=m_corner->getID();
}

/*!
  calculate & apply forces
*/
void BCorner2DInteraction::calcForces()
{
  // get number of adjacent edges
  int m_ne=m_corner->getNEdges();
  // transform anchor point to world coords
  Vec3 ap_global;
  if(m_ne==1){
  } else if (m_ne==2){
    ap_global=m_corner->getPos()+k1*m_corner->getEdgeNormal(1)+k2*m_corner->getEdgeNormal(2);
  }
  // get dist between anchor and particle
  const Vec3 D=ap_global-m_p->getPos();
  m_dist=sqrt(D*D);
  // calc force
  Vec3 force=D*m_k;
  Vec3 pos=m_p->getPos();
  // apply force
  m_p->applyForce(force,pos);
  if(m_ne==1){
  } else if (m_ne==2){
    Vec3 hf=force*(-0.5);
    m_corner->applyForceToEdge(1,hf);
    m_corner->applyForceToEdge(2,hf);
  }
}


/*!
  return if the interaction is broken, i.e. the distance between 
  particle and anchor point exceeds breaking distance, i.e. relative
  breaking distance x particle readius
*/
bool BCorner2DInteraction::broken()
{
  return (m_dist>m_break);
}

/*!
   Pack a BCorner2DInteraction into a TML packed message

  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::pack<BCorner2DInteraction>(const BCorner2DInteraction& I)
{
  append(I.m_k);
  append(I.m_dist);
  append(I.m_break);
  append(I.k1);
  append(I.k2);
  append(I.getCid());
  append(I.getPid());
}

/*!
  Unpack a BCorner2DInteraction from a TML packed message

  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::unpack<BCorner2DInteraction>(BCorner2DInteraction& I)
{
  I.m_k=pop_double();
  I.m_dist=pop_double();
  I.m_break=pop_double();
  I.k1=pop_double();
  I.k2=pop_double();
  I.m_cid=pop_int();
  I.m_pid=pop_int();
}
