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

#ifndef __FRAME_GEO_H
#define __FRAME_GEO_H

// -- STL includes --
#include <string>

using std::string;

/*!
  \author Feng Chen, Vince Boros, Michele Griffa
*/
void do_single_frame_geo(const string&,const string&,int,const string&);

/*!
  \author Feng Chen, Vince Boros, Michele Griffa
*/
void do_single_frame_geo_r(const string&,const string&,int,bool,const string&,bool,double,bool,const string&);

#endif // __FRAME_GEO_H
