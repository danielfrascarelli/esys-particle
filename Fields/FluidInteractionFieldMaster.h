/////////////////////////////////////////////////////////////
//                                                         //
// Copyright (c) 2003-2014 by The University of Queensland //
// Centre for Geoscience Computing                         //
// http://earth.uq.edu.au/centre-geoscience-computing      //
//                                                         //
// Primary Business: Brisbane, Queensland, Australia       //
// Licensed under the Open Software License version 3.0    //
// http://www.apache.org/licenses/LICENSE-2.0          //
//                                                         //
/////////////////////////////////////////////////////////////

#ifndef __FLUIDINTERACTIONFIELDMASTER_H
#define __FLUIDINTERACTIONFIELDMASTER_H

//--- project includes ---
#include "FieldMaster.h"
#include "Foundation/vec3.h"
#include "Foundation/quintuple.h"
#include "Foundation/triplet.h"
#include "field_const.h"

// --- TML includes ---
#include "tml/comm/comm.h"

//--- IO includes ---
#include <iostream>
#include <fstream>

//--- STL includes ---
#include <vector>
#include <string>
#include <map>

using std::vector;
using std::map;
using std::multimap;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;

class TML_Comm;

/*!
  \class ScalarFluidInteractionFieldMaster
  \brief Class for master part of a scalar field which is defined on all fluid interactions

  \author Qi Shao
  $Revision$
  $Date$
*/
class ScalarFluidInteractionFieldMaster : public AFieldMaster
{
 protected:
  vector<pair<Vec3,double> >  m_data_with_pos; // vector of <position,value> pairs
  vector<pair<int,double> >  m_data_with_id;
  vector<pair<pair<int,Vec3>,double> >  m_data_with_id_pos;
  vector<pair<esys::lsm::triplet<int,Vec3,double>,double> > m_data_with_particle;

  vector<double> m_sum_vec;
  vector<double> m_max_vec;

  virtual void writeAsRawWithPos();
  virtual void writeAsRawWithID();
  virtual void writeAsRawWithIDPos();
  virtual void writeAsRawWithParticle();
  virtual void writeAsSUM();
  virtual void writeAsMAX();
  virtual void writeAsRAW_SERIES();

  void collectFullwithPos();
  void collectFullwithID();
  void collectFullwithIDPos();
  void collectFullwithParticle();
  void collectSum();
  void collectMax();

 public:
  ScalarFluidInteractionFieldMaster(TML_Comm*,const string&,const string&,const string&,int,int,int);
  virtual ~ScalarFluidInteractionFieldMaster(){};

  void collect();
};


/*!
  \class VectorFluidInteractionFieldMaster
  \brief Class for master part of a vector field which is defined on all fluid interactions

  \author Qi Shao
  $Revision$
  $Date$
*/
class VectorFluidInteractionFieldMaster : public AFieldMaster
{
 protected:
  vector<pair<Vec3,Vec3> > m_data_with_pos; // vector of <position,value> pairs
  vector<pair<int,Vec3> >  m_data_with_id;
  vector<pair<pair<int,Vec3>,Vec3> > m_data_with_id_pos;
  vector<pair<esys::lsm::triplet<int,Vec3,double>,Vec3> > m_data_with_particle;

  vector<Vec3> m_sum_vec;

  virtual void writeAsRawWithPos();
  virtual void writeAsRawWithID();
  virtual void writeAsRawWithIDPos();
  virtual void writeAsRawWithParticle();
  virtual void writeAsSUM();
  virtual void writeAsRAW_SERIES();

  void collectFullwithPos();
  void collectFullwithID();
  void collectFullwithIDPos();
  void collectFullwithParticle();
  void collectSum();

 public:
  VectorFluidInteractionFieldMaster(TML_Comm*,const string&,const string&,const string&,int,int,int);
  virtual ~VectorFluidInteractionFieldMaster(){};

  void collect();
};
#endif //__FlUIDINTERACTIONFIELDMASTER_H
