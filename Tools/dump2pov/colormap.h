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

#ifndef __COLORMAP_H
#define __COLORMAP_H

#include "vec3.h"

class ColorMap
{
 protected:
  Vec3 c_min,c_max;
  double x_min,x_max;

 public:
  ColorMap(){};
  ColorMap(const Vec3&,const Vec3&,double,double);
  virtual Vec3 getColor(double) const;
};

#endif // __COLORMAP_H
