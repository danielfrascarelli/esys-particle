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

#ifndef __ECORNERINTERACTION_H
#define __ECORNERINTERACTION_H

// -- Project includes --
#include "Foundation/vec3.h"
#include "Geometry/Corner.h"
#include "Model/Particle.h"
#include "Model/ETriMeshIP.h"

/*!
  \class ECornerInteraction
  \brief unbonded elastic interaction between a Corner in a TriangleMesh and a particle

  \author Steffen Abe
  $Revision$
  $Date$
*/
class ECornerInteraction
{
 private:
  CParticle *m_p;
  Corner* m_corner;
  double m_k;//!< spring constant
   /*!
    flag showing if particle is in the inner area of the
    local particle array - needed for global force summation
  */
  bool m_inner_flag; 

 public:
  ECornerInteraction();
  ECornerInteraction(CParticle*,Corner*,ETriMeshIP,bool iflag=true);
  virtual ~ECornerInteraction();

  bool isInner(){return m_inner_flag;};
  virtual void calcForces();
};

#endif //__ECORNERINTERACTION_H
