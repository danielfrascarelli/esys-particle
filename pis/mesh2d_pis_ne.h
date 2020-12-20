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
#ifndef __MESH2D_PIS_NE_H
#define __MESH2D_PIS_NE_H

// --- project includes ---
#include "pis/mesh2d_pis.h"

/*!
  \class Mesh2D_PIS_NE
  \brief Class for parallel storage of interactions between a 2D 
  mesh and particles which doesn't require exchange of interactions across
  process boundaries

  \author Steffen Abe
  $Revision: $
  $Date: $
*/ 
template<class ParticleType,class IType>
class Mesh2D_PIS_NE : public Mesh2D_PIS<ParticleType>
{
 protected:
  typename IType::ParameterType m_param;

  set<pair<int,int> > m_edge_int_set; // for isIn
  set<pair<int,int> > m_corner_int_set; // for isIn
  vector<typename IType::EdgeIntType> m_edge_interactions;
  vector<typename IType::CornerIntType> m_corner_interactions;

 public:
  Mesh2D_PIS_NE(Mesh2D*,ParallelParticleArray<ParticleType>*,typename IType::ParameterType);
  ~Mesh2D_PIS_NE();

  virtual bool isIn(const vector<int>&);

  /**
   * Null op, time step size not required.
   */
  virtual void setTimeStepSize(double dt)
  {
  }

  virtual void calcForces();
  virtual bool update();
  virtual void exchange(){}; //!< do nothing
  virtual void rebuild(){}; //!< do nothing
  virtual void tryInsert(const vector<int>&){};//!< do nothing
};

#include "mesh2d_pis_ne.hpp"

#endif // __MESH2D_PIS_NE_H
