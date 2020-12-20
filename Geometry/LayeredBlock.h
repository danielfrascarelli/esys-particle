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

#ifndef __LAYEREDBLOCK_H
#define __LAYEREDBLOCK_H

//-- project includes --
#include "Geometry/RandomBlock.h"

//-- STL includes --
#include <set>

using std::set;

/*!
  \class CLayeredBlock2D
  \brief Class for the generation of a layered 2D random lattice in a rectangular area. 

  \author Steffen Abe
  $Revision$
  $Data:$
*/
class CLayeredBlock2D : public CRandomBlock2D
{
 private:
  set<double> LayerBoundaries;

 public:
  CLayeredBlock2D(double,double,double,double,double,double);
  virtual ~CLayeredBlock2D();

  void addLayerBoundary(double);
  virtual void generate(int,unsigned int);
};

#endif //__LAYEREDBLOCK_H
