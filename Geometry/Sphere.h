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

#ifndef __SPHERE2D_H
#define __SPHERE2D_H

//-- project includes --
#include "vec3.h"

/*! 
  \class Sphere2D
  \brief Just methods to calculate the position and dimension of a 2D sphere under given constraints (see the .cpp file).
  
  \author David Place, Steffen Abe
  $Revision$
  $Date$
*/ 

class Sphere {
public:
  // 2D fill-in
  static bool FillIn(Vec3 P1, Vec3 P2, Vec3 P3, double r1, double r2, double r3, Vec3 &M, double &r) ;
  static bool FillInWP(Vec3 P1, Vec3 P2, Vec3 O, Vec3 iD, double r1, double r2, Vec3 &M, double &r) ;
  static bool FillInWP(Vec3 iP1, Vec3 iO, Vec3 iD, double r1, double r, Vec3 &M, int wsol=1) ;
  //static bool FillInWFS(Vec3 P1, Vec3 P2, AGeneralSurface &FS, double r1, double r2, Vec3 &M, double &r) ;
  //static bool FillInWFS(Vec3 P1, AGeneralSurface &FS, double r1, double r, Vec3 &M, int wsol=1) ;
} ;

#endif // __SPHERE2D_H
