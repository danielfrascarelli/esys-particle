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

#include "dslice.h"

bool operator==(const DSlice& s1,const DSlice& s2)
{
  return (s1.m_start==s2.m_start &&
	  s1.m_length1==s2.m_length1 &&
	  s1.m_stride1==s2.m_stride1 &&
	  s1.m_length2==s2.m_length2 &&
	  s1.m_stride2==s2.m_stride2);
}

bool operator!=(const DSlice& s1,const DSlice& s2)
{
  return (s1.m_start!=s2.m_start ||
	  s1.m_length1!=s2.m_length1 ||
	  s1.m_stride1!=s2.m_stride1 ||
	  s1.m_length2!=s2.m_length2 ||
	  s1.m_stride2!=s2.m_stride2);
}

ostream& operator<< (ostream& ost,const DSlice& s)
{
  ost << "DSlice: " << endl;
  ost << "start: " << s.m_start << endl;
  ost << "length1, stride1 : " << s.m_length1 << ", " << s.m_stride1 << endl;
  ost << "length2, stride2 : " << s.m_length2 << ", " << s.m_stride2 << endl;

  return ost;
}
