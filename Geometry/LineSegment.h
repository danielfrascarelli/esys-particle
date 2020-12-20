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

#ifndef __LINE_SEGMENT_H
#define  __LINE_SEGMENT_H

//--- Project includes ---
#include "Foundation/vec3.h"
#include "Geometry/Line.h"

/*!
  \class LineSegment
  \brief Class representing a line segment for intersection/fitting calculation in 2D

  \author Steffen Abe
  $Date$
  $Revision$
*/
class LineSegment : public Line
{
 private:
  double m_len;

 public:
  LineSegment(const Vec3&,const Vec3&);
  virtual ~LineSegment(){}

  virtual double sep(const Vec3&);
  virtual bool intersect(const Vec3&,const Vec3&);
  Vec3 getP1() {return Pos;};
  Vec3 getP2() {return Pos+m_len*U;};
};

#endif // __LINE_SEGMENT_H
