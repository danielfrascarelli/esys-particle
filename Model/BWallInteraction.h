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

#ifndef __BWALLINTERACTION_H
#define __BWALLINTERACTION_H

#include "Model/WallInteraction.h"
#include "Model/Wall.h"

/*!
  \class CBondedWallInteraction
  \brief bonded elastic interaction between a particle and a wall

  \author Steffen Abe
  $Revision$
  $Date$  
*/

template <class T>
class CBondedWallInteraction : public AWallInteraction<T>
{
protected:
  double m_k;//!< spring constant
public:
  //  CBondedWallInteraction();
  CBondedWallInteraction(T*,CWall*,double,bool);
  virtual ~CBondedWallInteraction(){};

  virtual void calcForces();
  virtual Vec3 getForce();
  virtual void setPP(const vector<T*>){};
  virtual double getStiffness(){return m_k;};
};

#include "BWallInteraction.hpp"

#endif //__BWALLINTERACTION_H
