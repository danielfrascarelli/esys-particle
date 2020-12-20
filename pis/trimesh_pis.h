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

#ifndef __TRIMESH_PIS_H
#define __TRIMESH_PIS_H

// --- project includes ---
#include "pis/pi_storage.h"

// --- STL includes ---
#include <set>
#include <list>

/*!
  \brief Abstract base class for parallel storage of interactions between a triangle 
  mesh and particles
*/
template<class ParticleType> 
class TriMesh_PIS : public AParallelInteractionStorage
{
 protected:
  int m_update_timestamp;
  TriMesh* m_mesh;
 
 public:
  TriMesh_PIS(TriMesh*,ParallelParticleArray<ParticleType>*);
  virtual ~TriMesh_PIS();

  virtual void addExIG(AParallelInteractionStorage*);
  virtual AFieldSlave* generateNewScalarFieldSlave(TML_Comm*,const string&,int,int,int,int);
  virtual AFieldSlave* generateNewVectorFieldSlave(TML_Comm*,const string&,int,int,int,int);
};

#include "pis/trimesh_pis.hpp"

#endif //__TRIMESH_PIS_H
