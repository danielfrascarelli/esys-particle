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

#include "Line.h"

/*!
  protected empty constructor -> for use by derived classes (i.e. LineSegment) only
*/
Line::Line()
{}

Line::Line(const Vec3& D,const Vec3& P)
{
  N=D;
  U.X()=-N.Y();
  U.Y()=N.X();
  U.Z()=0.0;
  Pos=P;
}

double Line::sep(const Vec3& M)
{
  return fabs((M-Pos)*N);
}
