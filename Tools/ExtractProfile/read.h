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

#ifndef __READ_EXPROF_H
#define __READ_EXPROF_H

// -- STL includes --
#include <string>

using std::string;

void read_and_write_profile_r(const string&,const string&,double,double,int,bool,bool,int,int);
void read_and_write_profile_rel(const string&,const string&,const string&,double,double,int,bool,bool,int);
void read_and_write_disp_grid(const string&,const string&,double,double,double,double,double,double,double,bool,int,int);
void read_and_write_poros_grid(const string&,const string&,double,double,double,double,double,double,double);

#endif // __READ_EXPROF_H
