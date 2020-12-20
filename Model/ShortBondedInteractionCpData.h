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

#ifndef __SHORTBONDEDINTERACTIONCPDATA_H
#define __SHORTBONDEDINTERACTIONCPDATA_H

// --- Project includes ---
#include "Model/BondedInteractionCpData.h"

// --- IO includes ---
#include <iostream>

using std::istream;
using std::ostream;

class CShortBondedInteraction;

/*!
  Helper class for checkpointing ShortBondedInteraction data
*/ 
class  ShortBondedInteractionCpData : public BondedInteractionCpData
{
 private:
  double m_r0;
  
 public:
  ShortBondedInteractionCpData();
  ShortBondedInteractionCpData(int,int,int,double);
  ShortBondedInteractionCpData(const CShortBondedInteraction&);

  virtual ~ShortBondedInteractionCpData(){}
  
  virtual void saveCheckPointData(ostream&);
  virtual void loadCheckPointData(istream&);
};
#endif // __SHORTBONDEDINTERACTIONCPDATA_H
