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

#ifndef __SPHERE3D_H
#define __SPHERE3D_H

//-- project includes --
#include "Foundation/vec3.h"

/*! 
  \class Sphere3D
  \brief Just methods to calculate the position and dimension of a 2D sphere under given constraints (see the .cpp file).
  
  \author David Place, Steffen Abe
  $Revision$
  $Date$
*/ 

class Sphere3D {
 private:
  static double NearZero;
  
 public:
  static bool FillIn(const Vec3&,const Vec3&,const Vec3&,const Vec3&,double,double,double,double,Vec3&,double&);
  static bool FillInWP(const Vec3&,const Vec3&,const Vec3&,const Vec3&,const Vec3&,double,double,double,Vec3&,double&) ; 
};

#endif // __SPHERE3D_H
