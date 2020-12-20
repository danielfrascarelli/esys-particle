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

#include "InteractionParam.h"

//--- STL includes ---
#include <utility>
using std::pair;

AIParam::AIParam(const string& s)
{
  m_name=s;
}

AIParam::~AIParam()
{}

double AIParam::getParamByName(const string& s)
{
  return m_data[s];
}

void AIParam::addParameter(const string& s,double d)
{
  m_data.insert(pair<string,double>(s,d));
}
