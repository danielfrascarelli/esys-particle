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

//--- project includes ---
#include "Model/AdhesiveFriction.h"
#include "tml/message/packed_message_interface.h"

/*!
  Default constructor for CAdhesiveFriction interaction
  Zero all coefficients
*/
CAdhesiveFriction::CAdhesiveFriction()
{
  m_k=0.0;
  m_ks=0.0;
  m_r0=0.0;
  m_dt=0.0;
  m_r_cut=0.0;
  m_r_cut_h=0.0;
}

/*!
  Construct a CAdhesiveFriction interaction from 2 particle pointers and the parameters

  \param p1 pointer to the first particle
  \param p2 pointer to the second particle
  \param param the interaction parameters
*/
CAdhesiveFriction::CAdhesiveFriction(CParticle* p1,CParticle* p2,const CAdhesiveFrictionIGP& param):CFrictionInteraction(p1,p2)
{
  m_k=param.k;
  m_ks=param.k_s;
  m_r0=p1->getRad()+p2->getRad();
  m_dt=param.dt;
  m_r_cut=param.r_cut;
  m_r_cut_h=1.0+(param.r_cut-1.0)*0.5;  
}

/*!
  destruct a CAdehsiveFriction interaction, i.e.do nothing
*/
CAdhesiveFriction::~CAdhesiveFriction()
{}

void CAdhesiveFriction::calcForces()
{
  Vec3 pos;
  Vec3 force;
  
  // calculate distance
  Vec3 D=m_p1->getPos()-m_p2->getPos();
  double dist=D*D;
  // check if there is contact
  if(dist<(m_r0*m_r0)){ // contact -> calculate forces as for normal frictional interaction
    CFrictionInteraction::calcForces();
  } else if (dist<(m_r0*m_r0*m_r_cut_h*m_r_cut_h)){ // between eq and half cut -> increase
    //--- elastic force ---
    dist=sqrt(dist);
    force=D*(m_k*(dist-m_r0)/dist);
    pos=m_p2->getPos()+(m_p2->getRad()/m_r0)*D;
    m_Ffric=Vec3(0.0,0.0,0.0);
    m_normal_force=Vec3(0.0,0.0,0.0);
    // apply elastic force
    m_p2->applyForce(force,pos);
    m_p1->applyForce(-1.0*force,pos);
  } else if (dist<(m_r0*m_r0*m_r_cut*m_r_cut)){ // between half cut and cut -> decrease
    //--- elastic force ---
    dist=sqrt(dist);
    force=D*(m_k*((m_r0*m_r_cut)-dist)/dist);
    pos=m_p2->getPos()+(m_p2->getRad()/m_r0)*D;
    m_Ffric=Vec3(0.0,0.0,0.0);
    m_normal_force=Vec3(0.0,0.0,0.0);
    // apply elastic force
    m_p2->applyForce(force,pos);
    m_p1->applyForce(-1.0*force,pos);
  }
}

/*!
  Pack a CAdhesiveFriction into a TML packed message
 
  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::pack<CAdhesiveFriction>(const CAdhesiveFriction& I)
{
	append(I.m_k);
	append(I.m_r0);
	append(I.m_mu);
	append(I.m_ks);
	append(I.m_dt);
	append(I.m_r_cut);
	append(I.m_id[0]);
	append(I.m_id[1]);
}

/*!
  Unpack a CAdhesiveFriction from a TML packed message
 
  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::unpack<CAdhesiveFriction>(CAdhesiveFriction& I)
{
	I.m_k=pop_double();
	I.m_r0=pop_double();
	I.m_mu=pop_double();
	I.m_ks=pop_double();
	I.m_dt=pop_double();
	I.m_r_cut=pop_double();
	I.m_r_cut_h=1.0+(I.m_r_cut-1.0)*0.5; 
	I.m_id.erase(I.m_id.begin(),I.m_id.end());
	I.m_id.push_back(pop_int());
	I.m_id.push_back(pop_int());
}
