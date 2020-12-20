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

#ifndef __BEDGE2DINTERACTION_H
#define __BEDGE2DINTERACTION_H

// -- Project includes --
#include "Foundation/vec3.h"
#include "Geometry/Edge2D.h"
#include "Model/Particle.h"
#include "Model/BMesh2DIP.h"
#include "Model/BMesh2DInteractionCpData.h"

/*!
  \class BEdge2DInteraction
  \brief bonded elastic interaction between an edge in a 2d mesh and a particle

  \author Steffen Abe
  $Revision$
  $Date$
*/
class BEdge2DInteraction
{
 private:
  CParticle *m_p;
  Edge2D *m_ed;
  double m_k;
  double m_break;
  double m_dist;
  int m_eid;
  int m_pid;

  Vec3 m_ap; // anchor point in local coord. 
   /*!
    flag showing if particle is in the inner area of the
    local particle array - needed for global force summation
  */
  bool m_inner_flag; 

 public:  
  typedef BMesh2DIP  ParameterType;
  typedef BMesh2DInteractionCpData CheckPointable; //!< Used by PIS to save/load check-point data for objects of this type.

  BEdge2DInteraction();
  BEdge2DInteraction(CParticle*,Edge2D*,BMesh2DIP,bool iflag=true);
  virtual ~BEdge2DInteraction();

  bool isInner(){return m_inner_flag;};
  virtual void calcForces();
  bool broken();
  virtual Vec3 getPos()const {return m_p->getPos();}; // ??
  inline int getPid() const {return m_pid;};
  inline int getTid() const {return m_eid;};
  Vec3 getAP() const; //!< Return anchor point in global coordinates. Needed for snapshot data.
  virtual void setPP(CParticle* part_p){m_p=part_p;};
  virtual void setTP(Edge2D* tri_p){m_ed=tri_p;};

  friend class TML_PackedMessageInterface;
};
#endif //__BEDGE2DINTERACTION_H
