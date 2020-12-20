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

#include "Geometry/RectPatch.h"

// --- system includes ---
#include <cmath> // for fabs
using std::fabs;

/*!
  construct a axis aligned rectangular patch in the x-z plane

  \param xmin min. x-pos
  \param xmax max. x-pos
  \param zmin min. z-pos
  \param zmax max. z-pos
  \param z0 z-pos
  \param dz "roughness" parameter
*/
RectPatch::RectPatch(double xmin,double xmax,double zmin,double zmax,double y0,double dy)
{
  m_xmin=xmin;
  m_xmax=xmax;
  m_zmin=zmin;
  m_zmax=zmax;
  m_y0=y0;
  m_dy=dy;
}

/*!
  Get (perpendicular) distance from given point. 
  If projection of point onto plane is outside patch return -1

  \param P the point
*/
double RectPatch::sep(const Vec3& P)
{
  double res;

  // check if inside
  if((P.X()>=m_xmin) && (P.X()<=m_xmax) && (P.Z()>=m_zmin) && (P.Z()<=m_zmax)){
    res=fabs(P.Y()-m_y0)+m_dy;
  } else {
    res=-1;
  }

  return res;
}

/*!
  Get distance from given point to closest point of patch. 

  \param P the point
*/
double RectPatch::dist(const Vec3& P)
{
  return sep(P);
}

/*!
  check if line between 2 points intersects patch

  \param P1 1st point 
  \param P2 2nd point
*/
bool RectPatch::intersect(const Vec3& P1,const Vec3& P2)
{
  bool res;

  // check if patch is in z-range of line
  double dy1=P1.Y()-m_y0;
  double dy2=P2.Y()-m_y0;
  if(dy1*dy2<0.0){
    // get intersection point
    Vec3 P=P1+(dy1/(dy1+dy2))*(P2-P1);
    // check inside
    res=((P.X()>=m_xmin) && (P.X()<=m_xmax) && (P.Z()>=m_zmin) && (P.Z()<=m_zmax));
  } else { 
    res=false; 
  }

  return res;
}

/*!
  get the plane further away from the given point
*/
Plane3D RectPatch::getPlane(const Vec3& P)
{
  Vec3 normal=(P.Y()>m_y0) ? Vec3(0.0,1.0,0.0) :  Vec3(0.0,-1.0,0.0);
  Vec3 pos=(P.Y()>m_y0) ? Vec3(m_xmin,m_y0-m_dy,m_zmin) : Vec3(m_xmin,m_y0+m_dy,m_zmin) ;

  return Plane3D(normal,pos);
}
