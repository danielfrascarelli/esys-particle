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

#ifndef __TRIMESH_PIS_EB_H
#define __TRIMESH_PIS_EB_H

// --- STL includes ---
#include <vector>

// --- project includes ---
#include "tml/comm/cart_comm.h"
#include "pis/trimesh_pis.h"

/*!
  \brief Class for parallel storage of interactions between a triangle
  mesh and particles which does require exchange of interactions across
  process boundaries but where interactions are not dynamically formed
*/
template<class ParticleType,class IType>
class TriMesh_PIS_EB  : public TriMesh_PIS<ParticleType>
{
 private:
  static const int m_exchg_tag;
  void exchange_boundary(int,int);

 protected:
  typename IType::ParameterType m_param;

  TML_CartComm m_comm;
  set<pair<int,int> > m_tri_int_set; // for isIn, <TID,PID> pairs 
  list<typename IType::TriIntType> m_triangle_interactions;

 public:
  TriMesh_PIS_EB(TriMesh*,ParallelParticleArray<ParticleType>*,typename IType::ParameterType);
   
  virtual bool isIn(const vector<int>&);
  virtual void setTimeStepSize(double dt);
  virtual void calcForces();
  virtual bool update();
  virtual void exchange();
  virtual void rebuild();
  virtual void tryInsert(const typename IType::TriIntType&);
  virtual void tryInsert(const vector<int>&);

  virtual void saveSnapShotData(std::ostream&);

  void buildFromPPATagged(int,int);
  void buildFromPPAByGap(double);
};

#include "pis/trimesh_pis_eb.hpp"

#endif // __TRIMESH_PIS_EB_H
