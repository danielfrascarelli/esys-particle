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

#include "Geometry/PointCloud.h"

/*!
  construct an empty point cloud - do nothing
*/
PointCloud::PointCloud()
{}

/*!
  add a point 

  \param p the position of the point
*/
void PointCloud::addPoint(const Vec3& p)
{
  m_points.push_back(p);
}

/*!
  calculate center point
*/
Vec3 PointCloud::getCenter()
{
  Vec3 Center=Vec3(0.0,0.0,0.0);
  for(vector<Vec3>::iterator iter=m_points.begin();
      iter!=m_points.end();
      iter++){
    Center+=(*iter);
  }
  Center=Center/m_points.size();

  return Center;
}

/*!
  find a plane that best fits the could of points. Algorithm as described in Schneider & Eberly
  "Geometric Tools for Computer Graphics", pp 884/885, i.e. getting the eigenvectors of a matrix
  with 
  m_11=\sum(x_i-a)^2
  m_12=\sum(x_i-a)(y_i-b)
  m_13=\sum(x_i-a)(z_i-c)
  m_22=\sum(y_i-b)^2
  m_23=\sum(y_i-b)(z_i-c)
  m_33=\sum(z_i-c)^2
  where (a,b,c) is the center of the cloud
*/
Plane3D PointCloud::getFitPlane()
{
  // calculate center point
  Vec3 C=getCenter();

  // get matrix coefficients
  double m_11=0.0;
  double m_12=0.0;
  double m_13=0.0;
  double m_22=0.0;
  double m_23=0.0;
  double m_33=0.0;
  for(vector<Vec3>::iterator iter=m_points.begin();
      iter!=m_points.end();
      iter++){
    double dx=(iter->X()-C.X());
    double dy=(iter->Y()-C.Y());
    double dz=(iter->Z()-C.Z());
    m_11+=dx*dx;
    m_12+=dx*dy;
    m_13+=dx*dz;
    m_22+=dy*dy;
    m_23+=dy*dz;
    m_33+=dz*dz;
  }

  // setup matrix
  Matrix3 M;
  M(0,0)=m_11;
  M(0,1)=m_12;
  M(1,0)=m_12;
  M(0,2)=m_13;
  M(2,0)=m_13;
  M(1,1)=m_22;
  M(1,2)=m_23;
  M(2,1)=m_23;
  M(2,2)=m_33;
 
//   std::cout << "M: " << M << std::endl;
  // get eigenvectors
  Vec3 v1,v2,v3;
  double e1,e2,e3;

  M.eigen(v1,v2,v3,e1,e2,e3);
//   std::cout << "Eigenvalues : " << e1 << " , " << e2 << " , " << e3 << std::endl;
//   std::cout << "Eigenvectors: " << v1 << " , " << v2 << " , " << v3 << std::endl;

  Plane3D P(v1,C);
  
  return P;
}
