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

#ifndef __LINE_H
#define __LINE_H

//-- Project includes --
#include "Foundation/vec3.h"

/*!
  \class Line
  \brief Class representing a line

  \author Steffen Abe
*/
class Line
{
 protected:
  Vec3 Pos,U,N;
  Line();

 public:
  Line(const Vec3&,const Vec3&);
  virtual ~Line(){};

  Vec3 GetU() const {return U;};
  Vec3 GetO() const {return Pos;};
  Vec3 GetN() const {return N;};
  virtual double sep(const Vec3&);
};

#endif //__LINE_H
