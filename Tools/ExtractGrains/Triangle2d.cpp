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

#include "Triangle2d.h"
#include <iostream>

Triangle2D::Triangle2D(const Vec3& p0,const Vec3& p1,const Vec3& p2,int id):m_p0(p0),m_p1(p1),m_p2(p2),m_id(id)
{}

bool Triangle2D::isIn(const Vec3& P) const
{
  bool res;
  int icount=1;

  // get intersection with 0-1edge
  double dy0=P.Y()-m_p0.Y();
  double dy1=P.Y()-m_p1.Y();
  double dy2=P.Y()-m_p2.Y();
//   std::cout <<  "dy0,dy1,dy2 : " << dy0 << " , " << dy1 << " , " << dy2 << std::endl;
  if(dy0*dy1<0){
    // get x-coord. of intersection
    double x=m_p0.X()+dy0*((m_p1.X()-m_p0.X())/(m_p1.Y()-m_p0.Y()));
    //std::cout << "x0 : " << x << std::endl;
    if(P.X()>x) {
      icount=icount*2;
    } else {
      icount=icount*-2;
    } 
  } 
  // set intesection with 0-2 edge
  if(dy0*dy2<0){
    // get x-coord. of intersection
    double x=m_p0.X()+dy0*((m_p2.X()-m_p0.X())/(m_p2.Y()-m_p0.Y()));
    //    std::cout << "x1 : " << x << std::endl;
    if(P.X()>x) {
      icount=icount*2;
    } else {
      icount=icount*-2;
    } 
  } 
  // set intesection with 1-2 edge
  if(dy1*dy2<0){
    // get x-coord. of intersection
    double x=m_p1.X()+dy1*((m_p2.X()-m_p1.X())/(m_p2.Y()-m_p1.Y()));
    //    std::cout << "x2 : " << x << std::endl;
    if(P.X()>x) {
      icount=icount*2;
    } else {
      icount=icount*-2;
    } 
  } 

  // check result : if we have 2 intersections and the point is once left & once right , we get 1, -2, 2 -> -4
  if (icount==-4){
    res=true;
  } else {
    res=false;
  }

  return res;
}
