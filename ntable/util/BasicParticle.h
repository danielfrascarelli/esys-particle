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

#ifndef __BASICPARTICLE_H
#define __BASICPARTICLE_H

// -- project includes --
#include "vec3.h"

/*!
  \class CBasicParticle
  \brief Basic Particle class. 

  Contains only the "geometric part" of the particle, i.e. position and radius, no forces or such. 
  -- Modified for testing the new neighbortable ! ---

  \author Steffen Abe
  $Revision$
  $Date$
*/
class CBasicParticle
{
protected:
  Vec3 m_pos;
  double m_rad;
  int m_id;

public:
  CBasicParticle();
  CBasicParticle(int,const Vec3&,double);

  inline Vec3& getPPos() {return m_pos;};
  inline Vec3 getPos() const {return m_pos;};
  inline double getRad() const {return m_rad;};
  inline int getID() const {return m_id;};

  inline void moveBy(Vec3 v){m_pos+=v;};//!< move relative to current position
  inline void moveTo(Vec3 v){m_pos=v;}; //!< move absolute
  inline void setRad(double r){m_rad=r;};
};

#endif //__BASICPARTICLE_H
