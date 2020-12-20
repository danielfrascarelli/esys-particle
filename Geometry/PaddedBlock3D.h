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

#ifndef __PADDEDBLOCK3D_H
#define __PADDEDBLOCK3D_H

//-- project includes --
#include "Geometry/SplitBlock3D.h"

using std::vector;

/*!
  \class CPaddedBlock3D
  \brief Class for the generation of a 3D lattice with a random middle section
  in a rectangular area. 

  \author Steffen Abe
  $Revision$
  $Date$
*/  
class CPaddedBlock3D : public CSplitBlock3D
{
 protected:
  double m_pad_size;

  virtual Vec3 getAPoint();
  void generate_regular_padding();

 public:
  CPaddedBlock3D(double,double,double,double,double,double,double,double,double,double,int,bool circ_x=false);
  virtual ~CPaddedBlock3D(){};

  virtual void generate(int,unsigned int);  
};
#endif //__PADDEDBLOCK3D_H
