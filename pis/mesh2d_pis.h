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

#ifndef __MESH2D_PIS_H
#define __MESH2D_PIS_H

// --- project includes ---
#include "pi_storage.h"

// --- STL includes ---
#include <set>
#include <list>

using std::set;
using std::list;

/*!
  \class Mesh2D_PIS
  \brief Abstract base class for parallel storage of interactions between a 2D 
  mesh and particles

  \author Steffen Abe
  $Revision$
  $Date$
*/
template<class ParticleType> 
class Mesh2D_PIS : public AParallelInteractionStorage
{
 protected:
  int m_update_timestamp;
  Mesh2D* m_mesh;
 
 public:
  Mesh2D_PIS(Mesh2D*,ParallelParticleArray<ParticleType>*);
  virtual ~Mesh2D_PIS();

  virtual void addExIG(AParallelInteractionStorage*);
  virtual AFieldSlave* generateNewScalarFieldSlave(TML_Comm*,const string&,int,int,int,int);
  virtual AFieldSlave* generateNewVectorFieldSlave(TML_Comm*,const string&,int,int,int,int);
  
  virtual void saveCheckPointData(std::ostream&);
  virtual void loadCheckPointData(std::istream&);
};

#include "mesh2d_pis.hpp"

#endif //__MESH2D_PIS_H
