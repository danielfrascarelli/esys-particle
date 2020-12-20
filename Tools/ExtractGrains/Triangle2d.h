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

#ifndef __TRIANGLE2D_H
#define __TRIANGLE2D_H

// --- project includes ---
#include "Foundation/vec3.h"

class Triangle2D
{
 private:
  Vec3 m_p0,m_p1,m_p2;
  int m_id;

 public:
  Triangle2D(const Vec3&,const Vec3&,const Vec3&,int);

  int Id() const {return m_id;};
  bool isIn(const Vec3&) const;
};

#endif // __TRIANGLE2D_H
