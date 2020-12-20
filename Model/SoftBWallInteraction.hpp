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


#ifndef MODEL_SOFTBWALLINTERACTION_HPP
#define MODEL_SOFTBWALLINTERACTION_HPP

#include "SoftBWallInteraction.h"

template <class T>
CSoftBondedWallInteraction<T>::CSoftBondedWallInteraction(T* p,CWall* w,double normalK,double shearK,bool scaling,bool iflag):
  AWallInteraction<T>(p,w,iflag)
{
  double scale;

  if (scaling) {
    // scale stiffness to particle cross section
    if(CParticle::getDo2dCalculations()){ // 2D
//      scale=2.0*this->m_p->getRad();
      scale=1.0;
    } else { // 3D
//      scale=3.1415926536*this->m_p->getRad()*this->m_p->getRad();
      scale=3.1415926536*this->m_p->getRad();
    }
  }
  else {
     scale = 1.0;
  }

  m_normalK=normalK*scale;
  m_shearK=shearK*scale;
}

/*!
  calculate bonded elastic forces.
*/
template <class T>
void CSoftBondedWallInteraction<T>::calcForces()
{
  const Vec3 &n  = this->m_wall->getNormal();
  const Vec3 relDisplacement =
    (
      this->m_p->getTotalDisplacement()
      -
      this->m_wall->getTotalDisplacement()
    );
  
  const double normalDist = dot(relDisplacement, n)/(n.norm());
  const Vec3 normalForce = ((m_normalK*normalDist)/n.norm())*n;

  /*
   * Tangent displacement vector is relDisplacement point
   * projected onto plane.
   */
  const Vec3 tangentDisplacement =
    (relDisplacement - ((normalDist/(n.norm()))*n));
  const Vec3 tangentForce = m_shearK*tangentDisplacement;

  const Vec3 totalForce = normalForce + tangentForce;
  this->m_p->applyForce(-1.0*totalForce,this->m_p->getPos());
  if(this->m_inner_flag) this->m_wall->addForce(totalForce);
}

/*!
  calculate and return the bonded elastic force
*/
template <class T>
Vec3 CSoftBondedWallInteraction<T>::getForce()
{
  const Vec3 &n  = this->m_wall->getNormal();
  const Vec3 relDisplacement =
    (
      this->m_p->getTotalDisplacement()
      -
      this->m_wall->getTotalDisplacement()
    );

  const double normalDist = dot(relDisplacement, n)/(n.norm());
  const Vec3 normalForce = ((m_normalK*normalDist)/n.norm())*n;

  /*
   * Tangent displacement vector is relDisplacement point
   * projected onto plane.
   */
  const Vec3 tangentDisplacement =
    (relDisplacement - ((normalDist/(n.norm()))*n));
  const Vec3 tangentForce = m_shearK*tangentDisplacement;

  const Vec3 totalForce = normalForce + tangentForce;

  return totalForce;
}


#endif
