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

#ifndef __SOFTBWALLINTERACTION_H
#define __SOFTBWALLINTERACTION_H

#include "WallInteraction.h"
#include "Wall.h"

/*!
  \class CSoftBondedWallInteraction
  \brief bonded elastic interaction between a particle and a wall with different 
  spring constants in the normal and shear directions

  \author Steffen Abe
  $Revision$
  $Date$  
*/
template <class T>
class CSoftBondedWallInteraction : public AWallInteraction<T>
{
protected:
  double m_normalK,m_shearK; //!< directional spring constants
public:
  CSoftBondedWallInteraction();
  CSoftBondedWallInteraction(T*,CWall*,double,double,bool,bool);
  virtual ~CSoftBondedWallInteraction(){};

  virtual void calcForces();
  virtual Vec3 getForce();
  virtual void setPP(const vector<T*>){};
  virtual double getStiffness(){return m_normalK;};
};

#include "Model/SoftBWallInteraction.hpp"

#endif //__SOFTBWALLINTERACTION_H
