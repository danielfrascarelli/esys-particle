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


#ifndef MODEL_EWALLINTERACTION_HPP
#define MODEL_EWALLINTERACTION_HPP

/*!
  constructor for elastic interaction between particle & wall

  \param p pointer to the particle
  \param w pointer to the wall
  \param k spring constant
  \param iflag flag if the particle is in the inner part of the local NTable
*/
template <class T>
CElasticWallInteraction<T>::CElasticWallInteraction(T* p,CWall* w,double k,bool iflag):
  AWallInteraction<T>(p,w,iflag)
{
  double f=1.0;//this->m_p->getRad();; 
  // scale elastic param
    if(!CParticle::getDo2dCalculations()){
//      f*=this->m_p->getRad();
      f*=3.141592654*this->m_p->getRad();
    }
  m_k=f*k;
}

/*!
  calculate free elastic forces.
*/
template <class T>
void CElasticWallInteraction<T>::calcForces()
{

  double dist=(this->m_p->getPos()-this->m_wall->getOrigin())*this->m_wall->getNormal(); 
  //  console.XDebug() << "pos, rad, dist: " << this->m_p->getPos() << " " << this->m_p->getRad() << " " << dist << "\n";

  if(dist<this->m_p->getRad()){
    Vec3 force=m_k*(this->m_p->getRad()-dist)*this->m_wall->getNormal();
    Vec3 pos=this->m_p->getPos()-dist*this->m_wall->getNormal();

 /*   Vec3 disp = this->m_p->getPos()- this->m_p->getOldPos() ;
    Vec3 x_d = Vec3(1.0,0.0,0.0);
    Vec3 y_d = Vec3(0.0,1.0,0.0);
    Vec3 fy_vert = -m_k*disp.Y()*y_d;
    Vec3 fx_vert = -m_k*disp.X()*x_d;
    force = force+ fx_vert + fy_vert ;*/


// friction wall 
/*
     Vec3 x_d = Vec3(1.0,0.0,0.0);
     Vec3 y_d = Vec3(0.0,1.0,0.0);
     Vec3 fx_vert, fy_vert ;
     double miu = 1.0; 

    if ( this->m_p->getVel().X() >0.0)  fx_vert = - miu*force.norm()*x_d;
    else fx_vert =  miu*force.norm()*x_d;

    if ( this->m_p->getVel().Y() >0.0)  fy_vert = - miu*force.norm()*y_d;
    else fy_vert =  miu*force.norm()*y_d;
     
    force = force+ fx_vert + fy_vert ;
 */
    this->m_p->applyForce(force,pos);
    if(this->m_inner_flag) this->m_wall->addForce(-1.0*force);
  }
}

/*!
  calculate & return free elastic forces, don't apply them
*/
template <class T>
Vec3 CElasticWallInteraction<T>::getForce()
{
  Vec3 force=Vec3(0.0,0.0,0.0);
  double dist=(this->m_p->getPos()-this->m_wall->getOrigin())*this->m_wall->getNormal();
  if(dist<this->m_p->getRad()){
    force=m_k*(this->m_p->getRad()-dist)*this->m_wall->getNormal();
  }

  return -1.0*force;
}

/*!
  Get stiffness of the interaction. Returns spring constant (m_k) if
  interaction is in contact, 0.0 otherwise.
*/
template <class T>
double CElasticWallInteraction<T>::getStiffness()
{
  double res=0.0;
  double dist=(this->m_p->getPos()-this->m_wall->getOrigin())*this->m_wall->getNormal();
  if(dist<this->m_p->getRad()){
    res=m_k;
  }

  return res;
}


#endif
