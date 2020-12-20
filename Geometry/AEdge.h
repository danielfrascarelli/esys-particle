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

#ifndef __AEDGE_H
#define __AEDGE_H

//-- Project includes --
#include "Foundation/vec3.h"

/*!
  \class AEdge
  \brief abstract base class for edges in mesh (2D or 3D)

  \author Steffen Abe
  $Revision$
  $Date$
*/
class AEdge
{
 protected:
  Vec3 m_p0,m_p1;
  
 public:
  AEdge(const Vec3&,const Vec3&);
  virtual ~AEdge();

  double sep(const Vec3&) const;
  pair<bool,double> dist(const Vec3&) const ; // signed separation according to direction of the normal
  Vec3 getBoundingBoxMin() const; 
  Vec3 getBoundingBoxMax() const;   
};
#endif // __AEDGE_H
