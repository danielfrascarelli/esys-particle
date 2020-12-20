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

#ifndef __GEOCOLORMAP_H
#define __GEOCOLORMAP_H

#include "colormap.h"

// --- STL includes ---
#include <vector>
using std::vector;

class GeoColorMap : public ColorMap
{
 private:
  vector<double> m_bdry;

 public:
  GeoColorMap(const Vec3&,const Vec3&,double,double,int,double);
  virtual ~GeoColorMap()
  {
  }
  virtual Vec3 getColor(double) const;  
};

#endif // __GEOCOLORMAP_H
