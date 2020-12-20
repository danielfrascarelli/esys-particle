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

#ifndef __COLORMAP3_H
#define __COLORMAP3_H

#include "colormap.h"

// --- STL includes ---
#include <vector>
using std::vector;

class ColorMap3 : public ColorMap
{
 private:
  Vec3 m_c0,m_c1,m_c2;
  double m_x0,m_x1,m_x2;

 public:
  ColorMap3(const Vec3&,const Vec3&,const Vec3&,double,double,double);
  virtual ~ColorMap3(){}
  virtual Vec3 getColor(double) const;  
};

#endif // __COLORMAP3_H
