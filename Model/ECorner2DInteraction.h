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
#ifndef __ECORNER2DINTERACTION_H
#define __ECORNER2DINTERACTION_H

// -- Project includes --
#include "Foundation/vec3.h"
#include "Geometry/Corner2D.h"
#include "Model/Particle.h"
#include "Model/ETriMeshIP.h"

/*!
  \class ECornerInteraction
  \brief unbonded elastic interaction between a Corner in a TriangleMesh and a particle

  \author Steffen Abe
  $Revision: 430 $
  $Date: 2004-09-28 06:51:39 +0100 (Tue, 28 Sep 2004) $
*/
class ECorner2DInteraction
{
 private:
  CParticle *m_p;
  Corner2D* m_corner;
  double m_k;//!< spring constant
   /*!
    flag showing if particle is in the inner area of the
    local particle array - needed for global force summation
  */
  bool m_inner_flag; 

 public:
  ECorner2DInteraction();
  ECorner2DInteraction(CParticle*,Corner2D*,ETriMeshIP,bool iflag=true);
  ~ECorner2DInteraction();

  bool isInner(){return m_inner_flag;};
  virtual void calcForces();
};

#endif //__ECORNER2DINTERACTION_H
