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

#ifndef __WALLFIELDSLAVE_H
#define __WALLFIELDSLAVE_H

// --- project includes ---
#include "Fields/FieldSlave.h"

// --- STL includes ---
#include <vector>

using std::vector;

class CWall;

/*!
  \class AWallFieldSlave
  \brief Abstract base class for slave part of field defined on a Wall

  $Revision$
  $Date$
*/
class AWallFieldSlave : public AFieldSlave
{
 protected:
  vector<CWall*> m_wall;

 public:
  AWallFieldSlave(TML_Comm*);
  void addWall(CWall*);
};

#endif // __WALLFIELDSLAVE_H
