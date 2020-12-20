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
#include "colormap3.h"

ColorMap3::ColorMap3(const Vec3& c0,const Vec3& c1,const Vec3& c2,double x0 ,double x1 ,double x2)
{
  m_c0=c0;
  m_c1=c1;
  m_c2=c2;
  m_x0=x0;
  m_x1=x1;
  m_x2=x2;
}

Vec3 ColorMap3::getColor(double x) const  
{
  Vec3 res;

  if(x<m_x0){
    res=m_c0;
  } else if (x>m_x2){
    res=m_c2;
  } else {
    if(x<m_x1){
      res=m_c0+(m_c1-m_c0)*(x-m_x0)/(m_x1-m_x0);
    } else {
      res=m_c1+(m_c2-m_c1)*(x-m_x1)/(m_x2-m_x1);
    }
  }

  return res;
}
