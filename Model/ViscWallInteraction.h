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

#ifndef __VISCWALLINTERACTION_H
#define __VISCWALLINTERACTION_H

#include "Model/WallInteraction.h"
#include "Model/Wall.h"

/*!
  \class CViscWallInteraction
  \brief bonded elastic interaction between a particle and a wall

  \author Steffen Abe
  $Revision$
  $Date$  
*/
template <class T>
class CViscWallInteraction : public AWallInteraction<T>
{
protected:
  double m_nu; //!< viscosity

public:
  //  CViscWallInteraction();
  CViscWallInteraction(T*,CWall*,double,bool);
  virtual ~CViscWallInteraction(){};

  virtual void calcForces();
  virtual Vec3 getForce();
  virtual void setPP(const vector<T*>){};
};

#include "Model/ViscWallInteraction.hpp"

#endif //__VISCWALLINTERACTION_H
