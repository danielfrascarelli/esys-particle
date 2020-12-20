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

#ifndef __ROUGHPADDEDBLOCK3D_H
#define __ROUGHPADDEDBLOCK3D_H

//-- project includes --
#include "Geometry/PaddedBlock3D.h"
#include "Geometry/Plane3D.h"
#include "Geometry/RectPatch.h"

// --- STL includes ---
#include <vector>

/*!
  \class CRoughPaddedBlock3D
  \brief Class for the generation of a 3D lattice with a random middle section
  and random rough/smooth sections of the fault surface

  \author Steffen Abe
  $Revision:  $
  $Date: $
*/ 
class CRoughPaddedBlock3D : public CPaddedBlock3D
{
 protected:
  double m_rough_xres,m_rough_yres;
  double m_rough_depth; 
  double m_rough_prob;

  vector<RectPatch> m_fault; 

  virtual RectPatch getClosestPatch(const SimpleParticle&,double);
  virtual Plane3D getClosestPlane(const SimpleParticle&);
 
 public:
  CRoughPaddedBlock3D(double,double,double,double,double,double,double,double,double,double,bool circ_x=false);
  virtual ~CRoughPaddedBlock3D(){};

  void setRoughness(int,int,double,double);

  virtual bool checkAFit(const SimpleParticle&) ;
  //virtual void generate(int,unsigned int);  
  virtual void generate(int);
/*   virtual void tagSplit(int,int,double); */
};

#endif // __ROUGHPADDEDBLOCK3D_H
