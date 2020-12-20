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


#ifndef SIMPLEPARTICLE_H
#define SIMPLEPARTICLE_H

#include "Foundation/console.h"
#include "Foundation/vec3.h"
#include "Geometry/SimpleParticleData.h"

/**
 *
 */
class SimpleParticle : public esys::lsm::SimpleParticleData
{
public:
  static const SimpleParticle INVALID;

  inline SimpleParticle(const Vec3 &posn, double radius, int id=0, int tag=0);

  inline SimpleParticle(const SimpleParticle &p);

  inline SimpleParticle &operator=(const SimpleParticle &p);

  inline const Vec3 &getPos() const;
  inline void setPos(const Vec3 &pos);
  inline void moveTo(const Vec3 &v);
  inline void translateBy(const Vec3 &v);
  inline void moveBy(const Vec3 &v);
  inline void rotate(const Vec3 &rotation, const Vec3 &posn);
  inline double getRad() const;
  inline void setRad(double r);

  inline bool isValid() const;

  template <typename TmplVisitor>
  void visit(const TmplVisitor &visitor) const;

  template <typename TmplVisitor>
  void visit(TmplVisitor &visitor);
};

inline std::ostream& operator<<(std::ostream &oStream, const SimpleParticle &particle);

/*!
  \class ParticleComparer
  \brief Compares distance of 2 particles to a 3rd particle

  $Date$
  $Revision$
*/
class ParticleComparer
{
 private:
  const SimpleParticle *m_pParticle;
 public:
  /**
   * Construct a comparison object for distances to a given particle.
   *
   *@param particle
   */
  inline ParticleComparer(const SimpleParticle&);
  inline bool operator()(const SimpleParticle&, const SimpleParticle&) const;
  inline bool operator()(const SimpleParticle*, const SimpleParticle*) const; 
};

#include "Geometry/SimpleParticle.hpp"

#endif
