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

#ifndef __SPHEREINTERACTION_H
#define __SPHEREINTERACTION_H

#include "SphereBody.h"
#include "Interaction.h"
#include "Particle.h"

/*! 
  \class ASphereBodyInteraction
  \brief Abstract base for all interactions between a particle and a sphere body
  
  \author Dion Weatherley
  $Revision$
  $Date$  
*/
template <class T>
class ASphereBodyInteraction : public AInteraction
{
   protected:
      T *m_p;
      CSphereBody *m_sphere;
      /*!
        flag showing if particle is in the inner area of the
        local particle array - needed for global force summation
      */
      bool m_inner_flag;

   public:
      ASphereBodyInteraction(T*,CSphereBody*,bool iflag=true);
      virtual ~ASphereBodyInteraction(){};

      virtual bool hasTag(int,int) const;
      virtual Vec3 getPosFirst() const {return m_p->getPos();};

      inline bool isInner(){return m_inner_flag;};
      virtual void calcForces()=0;
      virtual double getStiffness(){return 0.0;};
};

#include "SphereBodyInteraction.hpp"

#endif // __SPHEREINTERACTION_H
