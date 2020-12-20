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

#ifndef __BTRIANGLEINTERACTION_H
#define __BTRIANGLEINTERACTION_H

// -- Project includes --
#include "Foundation/vec3.h"
#include "Geometry/Triangle.h"
#include "Model/Particle.h"
#include "Model/BTriMeshIP.h"
#include "Model/BTriMeshInteractionCpData.h"

/*!
  \class BTriangleInteraction
  \brief bonded elastic interaction between a Triangle and a particle

  \author Steffen Abe
  $Revision$
  $Date$
*/
class BTriangleInteraction
{
 private:
  CParticle *m_p;
  Triangle *m_t;
  double m_k;
  double m_break;
  double m_dist;
  int m_tid;
  int m_pid;

  Vec3 m_ap; // anchor point in local triangle coords.
   /*!
    flag showing if particle is in the inner area of the
    local particle array - needed for global force summation
  */
  bool m_inner_flag; 

 public:
  typedef BTriMeshIP ParameterType;
  typedef BTriMeshInteractionCpData CheckPointable;

  BTriangleInteraction();
  BTriangleInteraction(CParticle*,Triangle*,BTriMeshIP,bool iflag=true);
  virtual ~BTriangleInteraction();

  bool isInner(){return m_inner_flag;};
  virtual void calcForces();
  bool broken();
  virtual Vec3 getPos()const {return m_p->getPos();}; // ??
  inline int getPid() const {return m_pid;};
  inline int getTid() const {return m_tid;};
  Vec3 getAP() const;
  virtual void setPP(CParticle* part_p){m_p=part_p;};
  virtual void setTP(Triangle* tri_p){m_t=tri_p;};

  friend class TML_PackedMessageInterface;
};
#endif //__BTRIANGLEINTERACTION_H
