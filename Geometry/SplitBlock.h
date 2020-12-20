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

#ifndef __SPLITBLOCK_H
#define __SPLITBLOCK_H

//-- project includes --
#include "RandomBlock.h"

/*!
  \class CSplitBlock2D
  \brief Class for the generation of a split 2D random lattice in a rectangular area. 

  \author Steffen Abe
  $Revision$
  $Data:$
*/
class CSplitBlock2D : public CRandomBlock2D
{
 protected:
  double m_ysplit;

 public:
  CSplitBlock2D(double,double,double,double,double,double,double,bool circ_x=false);
  virtual ~CSplitBlock2D();

  virtual void generate(int,unsigned int);
  virtual void tagSplit(int,int,double);
};

#endif // __SPLITBLOCK_H
