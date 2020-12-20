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

#ifndef __WALLINTERACTION_H
#define __WALLINTERACTION_H

#include "Wall.h"
#include "Interaction.h"
#include "Particle.h"

/*! 
  \class AWallInteraction
  \brief Abstract base for all interactions between a particle and a wall
  
  \author Steffen Abe
  $Revision$
  $Date$  
*/
template <class T>
class AWallInteraction : public AInteraction
{
protected:
  T *m_p;
  CWall *m_wall;
  /*!
    flag showing if particle is in the inner area of the
    local particle array - needed for global force summation
  */
  bool m_inner_flag; 

public:

  AWallInteraction(T*,CWall*,bool iflag=true);

  virtual ~AWallInteraction(){};

  virtual bool hasTag(int,int) const;
  virtual Vec3 getPosFirst() const {return m_p->getPos();};

  inline bool isInner(){return m_inner_flag;};
  virtual void calcForces()=0;
  virtual double getStiffness(){return 0.0;};
};

#include "WallInteraction.hpp"

#endif //__WALLINTERACTION_H
