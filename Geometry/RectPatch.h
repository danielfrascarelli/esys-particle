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

#ifndef __RECT_PATCH_H
#define __RECT_PATCH_H

// --- project includes ---
#include "Foundation/vec3.h"
#include "Geometry/Plane3D.h"

/*!
 */
class RectPatch
{
 private:
  double m_xmin,m_xmax,m_zmin,m_zmax,m_y0,m_dy;

 public:
  RectPatch(double,double,double,double,double,double);

  double sep(const Vec3&);
  double dist(const Vec3&);
  bool intersect(const Vec3&,const Vec3&);
  Plane3D getPlane(const Vec3&);
  Vec3 getBasePoint() const {return Vec3(m_xmin,m_y0,m_zmin);};
};



#endif // __RECT_PATCH_H
