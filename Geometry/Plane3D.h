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

#ifndef __PLANE3D_H
#define __PLANE3D_H

//-- Project includes --
#include "Foundation/vec3.h"

/*!
  \class Plane3D
  \brief Class representing a Plane3D

  \author Steffen Abe
  $Revision$
  $Date$

*/
class Plane3D
{
 protected:
  Vec3 U,V ;
  bool force3D ;
  Vec3 Dir ;  
  Vec3 Pos ;
  
  void Create() ;
 public:
  Plane3D();
  Plane3D(const Vec3& iDir,const Vec3& iPos);
  Plane3D(const Vec3& iU,const Vec3& iV,const Vec3& iPos);
    
  virtual ~Plane3D() {}
    
  virtual double sep(const Vec3&) const;
  virtual double dist(const Vec3&) ; // signed separation according to Direction of the normal
  virtual Vec3 ToClosest(const Vec3& M) {return dist(M)*Dir; } ; // return the vector PM where P is the closest point to M on the surface.

  inline Vec3 GetU() const { return U; } ; // return U (for planes in a 2D space (ie. line)  U is in the space)
  inline Vec3 GetV() const { return V; } ; // V is null if this is a 2D space.
  inline const Vec3 &GetW() const { return Dir; } ; // The normal
  inline Vec3 getNormal() const { return Dir; } ;
  inline const Vec3 &GetO() const { return Pos; } ;
  inline Vec3 getPos() const { return Pos; } ;
} ;

namespace esys
{
  namespace lsm
  {
    typedef ::Plane3D Plane3D;
  }
}

#endif // __PLANE3D_H
