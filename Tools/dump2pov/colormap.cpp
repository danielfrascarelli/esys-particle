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

#include "colormap.h" 

ColorMap::ColorMap(const Vec3& c0,const Vec3& cm,double x0,double xm)
{
  c_min=c0;
  c_max=cm;
  x_min=x0;
  x_max=xm;
}

Vec3 ColorMap::getColor(double x) const
{
  Vec3 res;

  if(x<x_min){
    res=c_min;
  } else if (x>x_max){
    res=c_max;
  } else {
    res=c_min+(c_max-c_min)*(x-x_min)/(x_max-x_min);
  }

  return res;
}
