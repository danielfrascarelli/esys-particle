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

#ifndef __BTRIMESHINTERACTIONCPDATA_H
#define __BTRIMESHINTERACTIONCPDATA_H

// -- Project includes --
#include "Foundation/vec3.h"
#include "Parallel/CheckPointable.h"

class BTriangleInteraction; // forward decl.

/**
 * Helper class for checkpointing BTriangleInteraction data.
 */
class BTriMeshInteractionCpData : public esys::lsm::CheckPointable
{
 private:
  Vec3 m_ap;
  int m_tid;
  int m_pid;

 public:
  BTriMeshInteractionCpData();
  virtual ~BTriMeshInteractionCpData(){}

  BTriMeshInteractionCpData(const BTriangleInteraction&);
  void set(const BTriangleInteraction& );
  virtual void saveSnapShotData(std::ostream&);
  virtual void saveCheckPointData(std::ostream&);
  virtual void loadCheckPointData(std::istream&);
};

#endif //__BTRIMESHINTERACTIONCPDATA_H
