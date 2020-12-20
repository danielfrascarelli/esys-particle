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

#include "AEdge.h"

/*!
  construct Edge from corner coordinates. 

  \param v0 first corner
  \param v1 second corner
*/
AEdge::AEdge(const Vec3& v0,const Vec3& v1)
{
  m_p0=v0; 
  m_p1=v1;
}

AEdge::~AEdge()
{}

/*!
  Get min. corner of axis-aligned bounding box 
*/
Vec3 AEdge::getBoundingBoxMin() const
{
  return cmin(m_p0,m_p1);
}

/*!
  Get max. corner of axis-aligned bounding box 
*/
Vec3 AEdge::getBoundingBoxMax() const
{
  return cmax(m_p0,m_p1);
}

/*!
  get distance between point and closest point along edge (incl. corners)

  \param p the point
*/
double AEdge::sep(const Vec3& p) const
{
  double sep;

  Vec3 v=m_p1-m_p0;
  Vec3 vu=v.unit();
  double d=((p-m_p0)*vu);
  if((d>0.0)&(d*d<v.norm2())){ // closest point within edge
    sep=((p-m_p0)-d*vu).norm();
  } else { // closest point outside -> check corner distances
    double d1=(p-m_p0).norm();
    double d2=(p-m_p1).norm();
    sep=(d1<d2) ? d1 : d2;
  }

  return sep;
}

/*!
  Get perpendicular distance between point and edge.
  If the closest point on the supportung line is outside the 
  edge, the first component of the return value is "false", 
  otherwise "true"

  \param p the point
*/
pair<bool,double> AEdge::dist(const Vec3& p) const 
{
  bool is_in=false;
  double dist=0.0;

  Vec3 v=m_p1-m_p0;
  Vec3 vu=v.unit();
  double d=((p-m_p0)*vu);
  if((d>0.0)&(d*d<v.norm2())){
    dist=((p-m_p0)-d*vu).norm();
    is_in=true;
  }else{
    is_in=false;
  }
  return make_pair(is_in,dist);
}
