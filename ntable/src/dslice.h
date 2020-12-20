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

#ifndef __DSLICE_H
#define __DSLICE_H

//--- IO includes ---
#include <iostream>
using std::ostream;
using std::endl;

class DSlice
{
 private:
  int m_start;
  int m_length1,m_length2;
  int m_stride1,m_stride2;

 public:
  DSlice(int st,int l1,int s1,int l2,int s2):
    m_start(st),m_length1(l1),m_length2(l2),m_stride1(s1),m_stride2(s2)
    {}
  
  inline unsigned int size() const {return m_length1*m_length2;}
  inline int operator[](int idx) const
    {return m_start+(idx%m_length1)*m_stride1+(idx/m_length1)*m_stride2;}

  friend bool operator==(const DSlice&,const DSlice&);
  friend bool operator!=(const DSlice&,const DSlice&);
  
  friend ostream& operator<< (ostream&,const DSlice&);
};

#endif //__DSLICE_H
