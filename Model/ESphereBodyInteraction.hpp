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


#ifndef MODEL_ESPHEREBODYINTERACTION_HPP
#define MODEL_ESPHEREBODYINTERACTION_HPP

/*!
  constructor for elastic interaction between particle & sphere body

  \param p pointer to the particle
  \param w pointer to the sphere body
  \param k spring constant
  \param iflag flag if the particle is in the inner part of the local NTable
*/
template <class T>
CElasticSphereBodyInteraction<T>::CElasticSphereBodyInteraction(T* p,CSphereBody* w,double k,bool iflag):
  ASphereBodyInteraction<T>(p,w,iflag)
{
  double f=1.0;
  // scale elastic param
    if(!CParticle::getDo2dCalculations()){
      f*=3.141592654*this->m_p->getRad();
    }
  m_k=f*k;
}

/*!
  calculate free elastic forces.
*/
template <class T>
void CElasticSphereBodyInteraction<T>::calcForces()
{

  Vec3 pvec=this->m_p->getPos()-this->m_sphere->getCentre();
  Vec3 normalVec=pvec.unit();
  double dist=pvec.norm() - this->m_sphere->getRadius();

  if(dist<this->m_p->getRad()){
    Vec3 force=m_k*(this->m_p->getRad()-dist)*normalVec;
    Vec3 pos=this->m_p->getPos()-dist*normalVec;

    this->m_p->applyForce(force,pos);
    if(this->m_inner_flag) this->m_sphere->addForce(-1.0*force);
  }
}

/*!
  calculate & return free elastic forces, don't apply them
*/
template <class T>
Vec3 CElasticSphereBodyInteraction<T>::getForce()
{
  Vec3 force=Vec3(0.0,0.0,0.0);
  Vec3 pvec=this->m_p->getPos()-this->m_sphere->getCentre();
  Vec3 normalVec=pvec.unit();
  double dist=pvec.norm() - this->m_sphere->getRadius();

  if(dist<this->m_p->getRad()){
    force=m_k*(this->m_p->getRad()-dist)*normalVec;
  }

  return -1.0*force;
}

/*!
  Get stiffness of the interaction. Returns spring constant (m_k) if
  interaction is in contact, 0.0 otherwise.
*/
template <class T>
double CElasticSphereBodyInteraction<T>::getStiffness()
{
  double res=0.0;
  Vec3 pvec=this->m_p->getPos()-this->m_sphere->getCentre();
  double dist=pvec.norm() - this->m_sphere->getRadius();
  if(dist<this->m_p->getRad()){
    res=m_k;
  }

  return res;
}

#endif
