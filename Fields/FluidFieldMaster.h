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

#ifndef __FLUIDFIELDMASTER_H
#define __FLUIDFIELDMASTER_H

//--- project includes ---
#include "FieldMaster.h"
#include "Foundation/vec3.h"
#include "field_const.h"

// --- TML includes ---
#include "tml/comm/comm.h"

//--- IO includes ---
#include <iostream>
#include <fstream>

//--- STL inculdes ---
#include <string>
#include <map>


using std::map;
using std::multimap;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;


class TML_Comm;

bool sortOnZ(const pair<Vec3,double>, const pair<Vec3,double>); 
bool sortOnY(const pair<Vec3,double>, const pair<Vec3,double>); 
bool sortOnX(const pair<Vec3,double>, const pair<Vec3,double>);
vector<pair<Vec3,double> > sortVector(vector<pair<Vec3,double> >);

/*!
  \class ScalarFluidFieldMaster
  \brief Class for master part of a scalar field which is defined on all fluid cells

  \author Qi Shao
  $Revision$
  $Date$
*/
class ScalarFluidFieldMaster : public AFieldMaster
{
 protected:
  vector<pair<Vec3,double> > m_save_vector;
  vector<double> m_sum_vector;

  void collectFull();
  void collectSum();
  virtual void writeAsSUM();
  virtual void writeAsMAX();
  virtual void writeAsRAW_SERIES();
  virtual void writeAsRAW(); 
  virtual void writeAsVTI();

 public:
  ScalarFluidFieldMaster(TML_Comm*,const string&,const string&,const string&,int,int,int);
  virtual ~ScalarFluidFieldMaster(){};

  virtual void collect();

};

/*!
  \class VectorFluidFieldMaster
  \brief Class for master part of a vector field which is defined on all fluid cells

  \author Qi Shao
  $Revision$
  $Date$
*/
class VectorFluidFieldMaster : public AFieldMaster
{
 protected:
  vector<pair<Vec3,Vec3> > m_save_vector;

  virtual void writeAsSUM();
  virtual void writeAsMAX();
  virtual void writeAsRAW_SERIES();
  virtual void writeAsRAW();
  virtual void writeAsVTU();

 public:
  VectorFluidFieldMaster(TML_Comm*,const string&,const string&,const string&,int,int,int);

  virtual ~VectorFluidFieldMaster(){};

  void collect();
};


#endif //__FLUIDFIELDMASTER_H
