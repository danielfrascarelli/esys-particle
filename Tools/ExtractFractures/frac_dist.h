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

#ifndef _FRAC_DIST_H
#define _FRAC_DIST_H

//--- STL includes ---
#include <string>
#include <vector>

using std::string;
using std::vector;

//--- project includes ---
#include "fracframe.h" // for fdata

class FracDist
{
 private:

  vector<int> m_dist;
  double m_dmin,m_dmax,m_dx;
  int m_nbin,m_ntotal;

 public:

  FracDist(const vector<FracFrame::fdata>,double,double,int);
  void write(const string&);
};

#endif // _FRAC_DIST_H
