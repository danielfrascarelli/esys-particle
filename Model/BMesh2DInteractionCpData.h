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

#ifndef __BMESH2DINTERACTIONCPDATA_H
#define __BMESH2DINTERACTIONCPDATA_H

// --- STL includes ---
#include <iostream>
using std::istream;
using std::ostream;

// --- project includes ---
#include "Parallel/CheckPointable.h"
#include "Foundation/vec3.h"

class BEdge2DInteraction;

/*!
  \class BMesh2DInteractioncpData
  \brief helper class to checkpoint bonded mesh2d interactions

*/
class BMesh2DInteractionCpData : public esys::lsm::CheckPointable
{
 private:
  Vec3 m_ap;
  int m_pid;
  int m_tid;

 public:
  BMesh2DInteractionCpData();
  virtual ~BMesh2DInteractionCpData()
  {
  }
  
  BMesh2DInteractionCpData(const BEdge2DInteraction&);
  BMesh2DInteractionCpData(int,int);

  void set(const BEdge2DInteraction&);
  void set(int,int);
  int getPID();
  int getTID();
 
  virtual void saveSnapShotData(ostream&);
  virtual void saveCheckPointData(ostream&);
  virtual void loadCheckPointData(istream&);
};


#endif // __BMESH2DINTERACTIONCPDATA_H
