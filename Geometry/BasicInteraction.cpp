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


#include "Geometry/BasicInteraction.h"

using std::ostream;
using std::endl;

BasicInteraction::BasicInteraction(Id p1, Id p2, Tag tag)
  : m_p1(p1),
    m_p2(p2),
    m_tag(tag)
{
}

ostream& operator<<(ostream& ost,const BasicInteraction& BI)
{
  ost << BI.m_p1 << " " << BI.m_p2 << " " << BI.m_tag;
  return ost;
}

bool BILess::operator()(const BasicInteraction& B1,const BasicInteraction& B2) 
{
  bool res;
  
  res=((B1.m_p1<B2.m_p1)||((B1.m_p1==B2.m_p1)&&(B1.m_p2<B2.m_p2)));
    
  return res;
}
