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

#include "Foundation/vec3.h"
#include "Foundation/vec3.hpp"

const Vec3 Vec3::ZERO = Vec3(0.0, 0.0, 0.0);

bool Vec3::operator<(const Vec3& rhs) const
{
  bool res;

  if(data[0]!=rhs.data[0]) {
    res=data[0]<rhs.data[0];
  } else if(data[1]!=rhs.data[1]){
    res=data[1]<rhs.data[1];
  } else if(data[2]!=rhs.data[2]){
    res=data[2]<rhs.data[2];
  } else {
    res=false;
  }

  return res;
}
