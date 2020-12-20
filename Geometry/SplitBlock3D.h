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

#ifndef __SPLITBLOCK3D_H
#define __SPLITBLOCK3D_H

//-- project includes --
#include "Geometry/RandomBlock3D.h"


/*!
  \class CSplitBlock3D
  \brief Class for the generation of a split 3D random lattice in a rectangular area. 

  \author Steffen Abe
  $Revision$
  $Data:$
*/
class CSplitBlock3D : public CRandomBlock3D
{
 protected:
  double m_ysplit;
  int m_dir;

 public:
  CSplitBlock3D(double,double,double,double,double,double,double,double,double,int,bool circ_x=false,bool rough=false);
  virtual ~CSplitBlock3D();

  virtual void generate(int,unsigned int);
  virtual void tagSplit(int,int,double);
};

#endif // __SPLITBLOCK3D_H
