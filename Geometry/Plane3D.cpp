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

#include "Geometry/Plane3D.h"

/*!
  Construct plane from pos. and normal vector

  \param iPos position
  \param iDir normal vector
*/ 
Plane3D::Plane3D(const Vec3& iDir,const Vec3& iPos) 
{
  Dir  = iDir.unit() ;
  Pos = iPos ; 
  Create() ;
}

/*!
  Construct plane from pos. and 2 vectors spanning the plane

  \param iPos position
  \param iU  
  \param iV
*/
Plane3D::Plane3D(const Vec3& iU,const Vec3& iV,const Vec3& iPos) 
{
  U=iU.unit();
  V=iV.unit(); 
  Vec3 Dir2=cross(U,V);
  Pos=iPos ;
  if ( (U*V) != 0.0 ) {
    V=U - ((U*U)/(U*V))*V ;
    V=V.unit();
  } else {
    V = iV.unit();
  }
  Dir= cross(U,V);
  if (Dir*Dir2 <0) {
    V *= -1.0 ;
    Dir *= -1.0 ;
  }
  Pos = iPos; 
}

/*!
  "empty" default constructor
*/
Plane3D::Plane3D() 
{
  U = Vec3(0.0,0.0,0.0) ;
  V = Vec3(0.0,0.0,0.0) ;
}

/*!
  setup spanning vectors from pos & normal
*/
void Plane3D::Create() 
{
  U = Vec3(0.0,0.0,0.0) ;
  V = Vec3(0.0,0.0,0.0) ;
  Vec3 X ; 
  X = Vec3(1.0,0.0,0.0) ;
  if ((cross(X,Dir)).norm2() == 0.0) X = Vec3(0.0,1.0,0.0) ;
  if ((cross(X,Dir)).norm2() == 0.0) X = Vec3(0.0,0.0,1.0) ;
  if ((cross(X,Dir)).norm2() != 0.0) {
    U = X - ((X*Dir)/(Dir*Dir))*Dir ;
    U = U.unit() ;
    V = cross(Dir,U) ;
  }
} 

/*!
  distance of a point from the plane

  \param M the point
*/
double Plane3D::sep(const Vec3& M) const
{
  return fabs((M-Pos)*Dir);
} 

/*!
  signed separation according to Direction of the normal

  \param M the point
*/
double Plane3D::dist(const Vec3& M) 
{
  return (M-Pos)*Dir ;
} 
