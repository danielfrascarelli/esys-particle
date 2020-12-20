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

#ifndef __ASPHEREBODYINTERACTIONGROUP_H
#define __ASPHEREBODYINTERACTIONGROUP_H

//--- project includes ---
#include "Model/SphereBody.h"
#include "Model/InteractionGroup.h"

//--- IO includes ---
#include <iostream>

//--- TML includes ---
#include "tml/comm/comm.h"

/*!
  \brief Abstract Base class for a group of interactions between particles and a sphere body
*/
template<class T>
class ASphereBodyInteractionGroup : public AInteractionGroup<T>
{
protected:
  CSphereBody* m_sphere; //!< the sphere body
  TML_Comm* m_comm; //!< MPI communicator
  int m_inner_count;

public:
  ASphereBodyInteractionGroup(TML_Comm* comm)
    :m_sphere(NULL),
     m_comm(comm),
     m_inner_count(0)
  {
  }

  virtual ~ASphereBodyInteractionGroup()
  {
  }

  /**
   * Null op, current sphere body interactions don't require time step size.
   */
  virtual void setTimeStepSize(double dt)
  {
    // do nothing, time step size not required in sphere body interactions.
  }

  virtual void calcForces()=0;
  
  virtual void applyForce(const Vec3&){
    std::cerr
      << "calling unimplemented function ASphereBodyInteractionGroup::applyForce"
      << std::endl;
  }
  virtual void setVelocity(const Vec3&){
    std::cerr
      << "calling unimplemented function ASphereBodyInteractionGroup::setVelocity"
      << std::endl;
  }
  inline double getDisplacement(){return m_sphere->getDisplacement();};
  inline void resetDisplacement(){m_sphere->resetDisplacement();};
  inline void moveSphereBodyBy(const Vec3& mv){m_sphere->moveBy(mv);};
  inline void zeroForce(){m_sphere->zeroForce();};
};


#endif // __ASPHEREBODYINTERACTIONGROUP_H
