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

#include "Model/IGParam.h"

AIGParam::AIGParam(const std::string &name) : m_name(name)
{
}

AIGParam::~AIGParam()
{
}

void AIGParam::packInto(CVarMPIBuffer* B) const
{
  B->append(m_name.c_str());
}

void AIGParam::setName(const std::string &name)
{
  m_name = name;
}
