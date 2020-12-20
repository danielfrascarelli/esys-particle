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

// --- project includes ---
#include "DataParticle.h"
#include "ntable.h"

/*!
  \class DataExtractor
  \brief class for the extraction of data from snapshots
*/
class DataExtractor 
{
 private:
  NeighborTable<DataParticle> m_data;

 public:
  DataExtractor(int,int,int,double,const Vec3&,const Vec3&); //modified (fluid contents)

  // I/O
  void read(const string&);
  void writeTensorDataVtk(const string&,const string&);
  void writeScalarDataVtk(const string&,const string&);
  
  // data extraction
  void StrainToTensorData(double);
  void MaxShearToScalarData();
};
