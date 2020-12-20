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

#ifndef __EEDGEINTERACTION_H
#define __EEDGEINTERACTION_H

// -- Project includes --
#include "Foundation/vec3.h"
#include "Geometry/Edge.h"
#include "Model/Particle.h"
#include "Model/ETriMeshIP.h"
/*!
  \class EEdgeInteraction
  \brief unbonded elastic interaction between a Edge in a TriangleMesh and a particle

  \author Steffen Abe
  $Revision$
  $Date$
*/
class EEdgeInteraction
{
 private:
  CParticle *m_p;
  Edge *m_edge;
  double m_k; //!< spring constant
   /*!
    flag showing if particle is in the inner area of the
    local particle array - needed for global force summation
  */
  bool m_inner_flag; 

 public:
  EEdgeInteraction();
  EEdgeInteraction(CParticle*,Edge*,ETriMeshIP,bool iflag=true);
  virtual ~EEdgeInteraction();

  bool isInner(){return m_inner_flag;};
  virtual void calcForces();
};
#endif // __EEDGEINTERACTION_H
