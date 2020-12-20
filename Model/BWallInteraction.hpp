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

#ifndef MODEL_BWALLINTERACTION_HPP
#define MODEL_BWALLINTERACTION_HPP

template <class T>
CBondedWallInteraction<T>::CBondedWallInteraction(T* p,CWall* w,double k,bool iflag):
  AWallInteraction<T>(p,w,iflag)
{
  // scale stiffness to particle cross section
  double scale;
  if(CParticle::getDo2dCalculations()){ // 2D
//    scale=2.0*this->m_p->getRad();
    scale=1.0;
  } else { // 3D
//    scale=3.1415926536*this->m_p->getRad()*this->m_p->getRad();
    scale=3.1415926536*this->m_p->getRad();
  }

  m_k=k*scale;
}

/*!
  calculate bonded elastic forces.
*/

template <class T>
void CBondedWallInteraction<T>::calcForces()
{
  Vec3 D=(this->m_p->getTotalDisplacement()-this->m_wall->getTotalDisplacement());
  //double dist=sqrt(D*D);
 
  Vec3 force=D*m_k;
  Vec3 pos=this->m_p->getPos();

  this->m_p->applyForce(-1.0*force,pos);
  if(this->m_inner_flag) this->m_wall->addForce(force);
}

/*!
  calculate and return the bonded elastic force
*/
template <class T>
Vec3 CBondedWallInteraction<T>::getForce()
{
  Vec3 D=(this->m_p->getTotalDisplacement()-this->m_wall->getTotalDisplacement());
  //const double dist=sqrt(D*D);
 
  return D*m_k;
}
 
#endif
