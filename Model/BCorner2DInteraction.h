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

#ifndef __BCORNER2DINTERACTION_H
#define __BCORNER2DINTERACTION_H

// --- project includes ---
#include "Geometry/Corner2D.h"
#include "Model/Particle.h"
#include "Model/BMesh2DIP.h"

/*!
  \class BCorner2DInteraction
  \brief bonded elastic interaction between Corner2D in a 2d mesh and a particle

  \author Steffen Abe
  $Revision$
  $Date$
*/
class BCorner2DInteraction
{
 private:
  CParticle *m_p;
  Corner2D* m_corner;
  double m_k;//!< spring constant
  double m_break;
  double m_dist;
  double k1,k2; //!< coefficients for calculating the anchor point from the normals of the adjacent edges
  int b_me;
  int m_cid;
  int m_pid;

  bool m_inner_flag;

 public:
  typedef BMesh2DIP  ParameterType;

  BCorner2DInteraction();
  BCorner2DInteraction(CParticle*,Corner2D*,BMesh2DIP,bool iflag=true);
  virtual ~BCorner2DInteraction(){};

  bool isInner(){return m_inner_flag;};
  virtual void calcForces();
  bool broken();
  virtual Vec3 getPos()const {return m_p->getPos();}; // ??
  inline int getPid() const {return m_pid;};
  inline int getCid() const {return m_cid;};
  virtual void setPP(CParticle* part_p){m_p=part_p;};
  virtual void setCP(Corner2D* corner_p){m_corner=corner_p;};

  friend class TML_PackedMessageInterface;
};

#endif // __BCORNER2DINTERACTION_H
