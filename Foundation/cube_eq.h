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

#ifndef __CUBE_EQ_H
#define __CUBE_EQ_H


// --- STL includes ---
#include <set>

using std::set;

/*!
  \class CubicEquation 
  \brief A class for a cubic equation. Used for eigenvalue calculation on 3D matrices

  \author Steffen Abe
  $Revision$
  $Date$
*/
class CubicEquation
{
 private:
  double m_a,m_b,m_c;

  double bisect(double,double,double);
  double f(double);

 public:
  CubicEquation(double,double,double);

  set<double> getRealRoots(double);
};

#endif //__CUBE_EQ_H
