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

#ifndef __AWALLINTERACTIONGROUP_H
#define __AWALLINTERACTIONGROUP_H

//--- project includes ---
#include "Model/Wall.h"
#include "Model/InteractionGroup.h"

//--- IO includes ---
#include <iostream>

//--- TML includes ---
#include "tml/comm/comm.h"

/*!
  \brief Abstract Base class for a group of interactions between particles and a wall
*/
template<class T>
class AWallInteractionGroup : public AInteractionGroup<T>
{
protected:
  CWall* m_wall; //!< the wall
  TML_Comm* m_comm; //!< MPI communicator
  int m_inner_count;

public:
  AWallInteractionGroup(TML_Comm* comm)
    :m_wall(NULL),
     m_comm(comm),
     m_inner_count(0)
  {
  }

  virtual ~AWallInteractionGroup()
  {
  }

  /**
   * Null op, current wall interactions don't require time step size.
   */
  virtual void setTimeStepSize(double dt)
  {
    // do nothing, time step size not required in wall interactions.
  }

  virtual void calcForces()=0;
  
  virtual void applyForce(const Vec3&){
    std::cerr
      << "calling unimplemented function AWallInteractionGroup::applyForce"
      << std::endl;
  }
  virtual void setVelocity(const Vec3&){
    std::cerr
      << "calling unimplemented function AWallInteractionGroup::setVelocity"
      << std::endl;
  }
  inline double getDisplacement(){return m_wall->getDisplacement();};
  inline void resetDisplacement(){m_wall->resetDisplacement();};
  inline void moveWallBy(const Vec3& mv){m_wall->moveBy(mv);};
  inline void setWallNormal(const Vec3& wn){m_wall->setNormal(wn);};
  inline void zeroForce(){m_wall->zeroForce();};
};


#endif // __AWALLINTERACTIONGROUP_H
