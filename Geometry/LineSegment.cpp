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

#include "LineSegment.h"

/*!
  constructor

  \param P0 1st end point
  \param P1 2nd end point

  \warning doesn't check P0!=P1
*/ 
LineSegment::LineSegment(const Vec3& P0,const Vec3& P1)
{
  Pos=P0;
  U=(P1-P0).unit();
  N.X()=U.Y();
  N.Y()=-U.X();
  N.Z()=0.0;
  m_len=(P1-P0).norm();
}

/*!
  distance between a point and the line segment

  \param P the position of the point
*/
double LineSegment::sep(const Vec3& P)
{
  double res;

  double du=(P-Pos)*U;
  if((0<=du) && (du <=m_len)){// nearest point inside segment
    res=fabs((P-Pos)*N); 
  } else { // nearest point  outside -> get distance to closest endpoint
    double d1=(P-Pos).norm();
    double d2=(P-(Pos+m_len*U)).norm();
    res = (d1<d2) ? d1 : d2;
  }

  return res;
}

/*!
  returns if the connecting line between two points intersects the line segment

  \param P1 1st point
  \param P2 2nd point
*/
bool LineSegment::intersect(const Vec3& P1,const Vec3& P2)
{
  bool res=false;

  Vec3 U2=(P2-P1).unit();

  // setup params
  double x1=U.X();
  double y1=U.Y();
  double x2=U2.X();
  double y2=U2.Y();
  double c=P1.X()-Pos.X();
  double d=P1.Y()-Pos.Y();

  // solve EQS
  double s=y2*x1-x2*y1;

  if(s!=0.0){ // if s==0 -> parallel segments -> no intersection
    double b=(c*y1-d*x1)/s;
    double a=(c*y2-d*x2)/s;
    res=((0<=a) && (a<=m_len) && (0<=b) && (b<=(P2-P1).norm()));
  } 

  return res;
}
