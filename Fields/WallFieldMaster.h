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

#ifndef __WALLFIELDMASTER_H
#define __WALLFIELDMASTER_H

// --- project includes ---
#include "FieldMaster.h"
#include "vec3.h"

// --- STL includes ---
#include <vector>
#include <map>

using std::vector;
using std::map;

class TML_Comm;

/*!
  \class VectorWallFieldMaster
  \brief Master part of vector field on walls

  $Revision$
  $Date$
*/
class VectorWallFieldMaster : public AFieldMaster
{
 protected:
  virtual void writeAsRAW_SERIES();
  virtual void writeAsSILO();
  map<int,Vec3> m_data;
  int m_sum_flag; // decides if field is sum of all slaves (force...) or not (pos...)

 public:
  VectorWallFieldMaster(TML_Comm*,const string&,const string&,vector<string>,const string&,int,int,int);
  virtual ~VectorWallFieldMaster(){};

  virtual void collect();
};
#endif // __WALLFIELDMASTER_H
