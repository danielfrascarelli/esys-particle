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
#ifndef __EEDGE2DINTERACTION_H
#define __EEDGE2DINTERACTION_H

// -- Project includes --
#include "Foundation/vec3.h"
#include "Geometry/Edge2D.h"
#include "Model/Particle.h"
#include "Model/ETriMeshIP.h"

/*!
  \class EEdgeInteraction
  \brief unbonded elastic interaction between a Edge in a TriangleMesh and a particle

  \author Steffen Abe
  $Revision: $
  $Date: $
*/
class EEdge2DInteraction
{
 private:
  CParticle *m_p;
  Edge2D *m_edge;
  double m_k; //!< spring constant
   /*!
    flag showing if particle is in the inner area of the
    local particle array - needed for global force summation
  */
  bool m_inner_flag; 

 public:
  EEdge2DInteraction();
  EEdge2DInteraction(CParticle*,Edge2D*,ETriMeshIP,bool iflag=true);
  virtual ~EEdge2DInteraction();

  bool isInner(){return m_inner_flag;};
  virtual void calcForces();
};
#endif // __EEDGE2DINTERACTION_H
