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
#include "Foundation/vec3.h"

// --- STL includes ---
#include <map>
#include <utility>

using std::map;
using std::pair;
using std::make_pair;

namespace esys
{
  namespace lsm
  {
    class SimpleParticleData;
  }
}

/*!
  \brief Basic Particle class. Contains only the "geometric part" of the
  particle, i.e. position and radius, no forces or such.

  \author Steffen Abe
  $Revision$
  $Date$

*/
class CBasicParticle
{
protected:
  Vec3 m_pos; //!< position
  double m_rad; //!< radius
  int m_global_id;
  int m_tag;

public:
  static const CBasicParticle INVALID;

  CBasicParticle();
  CBasicParticle(const Vec3 &pos, double radius, int id=-1, int tag=-1);
  CBasicParticle(const esys::lsm::SimpleParticleData &data);

  inline virtual ~CBasicParticle(){}

  inline Vec3 & getPPos() {return m_pos;}
  inline Vec3 getPos() const {return m_pos;}
  inline void setPos(const Vec3 &pos) {m_pos = pos;}
  inline double getRad() const {return m_rad;}
  inline int getID() const {return m_global_id;}
  inline void setID(int id) {m_global_id = id;}

  inline void moveBy(Vec3 v){m_pos+=v;} //!< move relative to current position
  inline void moveTo(Vec3 v){m_pos=v;} //!< move absolute
  inline void setRad(double r){m_rad=r;}

  //! particle tag handling
  inline void setTag(int t){m_tag=t;}
  inline int getTag() const {return m_tag;}
  inline bool isValid() const {return (getID() >= 0);}
};
ostream& operator<<(ostream&,const CBasicParticle&);

#endif //__BASICPARTICLE_H
