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

#ifndef __FAULTEDBLOCK2D_H
#define __FAULTEDBLOCK2D_H

//--- STL includes ---
#include <vector>
#include <utility>

using std::vector;
using std::pair;

//-- project includes --
#include "Geometry/RandomBlock.h"
#include "Geometry/LineSegment.h"

/*!
  \class FaultedBlock2D
  \brief class for the generation of a 2D random block with a fault consisting of line segments

  \author Steffen Abe
  $Date$
  $Revision$
*/
class FaultedBlock2D : public CRandomBlock2D
{
 protected:
  vector<pair<double,LineSegment> > m_fault;
  vector<LineSegment> m_f2;
  virtual Line *getClosestPlane(const SimpleParticle&);
  double m_pad_size;

  virtual Vec3 getAPoint();

 public:
  FaultedBlock2D(double,double,double,double,double,double,double,bool circ_x=false);
  virtual ~FaultedBlock2D();

  void addSegment(const Vec3&,const Vec3&,double);
  virtual bool checkAFit(const SimpleParticle&) ;
  virtual void generate(int,unsigned int);
  virtual void tagSplit(int,int,double);
};

#endif // __FAULTEDBLOCK2D_H
