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

#ifndef __SLIP2VTK2D_H
#define __SLIP2VTK2D_H

// --- STL includes ---
#include <string>

using std::string;

// 1D displacement data to 2D x-x_0(x,t) 
void slip2vtk2d(const string&,const string&,const string&,int,int,int,double,double,double,int);
// 1D displacement data to 2D dx(x,t) 
void slip2vtk2d_rate(const string&,const string&,const string&,int,int,int,double,double,double,int);
// 1D displacement data to 2D dx(x,t) using RAW_SERIES and initial pos. file
void slip2vtk2d_rate_rs(const string&, const string&, const  string&, const string&,  const string&, int,int,int,double,double,double,int);
// 1D displacement data to 2D dx(x,t) using RAW_SERIES and initial pos. file - output format raw
void slip2raw2d_rs(const string&, const string&, const  string&, const string&,  const string&, int,int,int,double,double,int,int,int);
// 1D displacement data to 2D dx(x,t) using series of RAW2  displacement files- output format raw
void slip2raw2d(const string&, const string&, const string&, int,int,int,double,double,int);
// 1D displacement data to 2D v(x,t) using RAW_SERIES and initial pos. file - output format raw
void slip2raw2d_rate_rs(const string&, const string&, const  string&, const string&,  const string&, int,int,int,double,double,int);
// 1D displacement data 1D dx(t) (x=x_0)
void slip_x2slip_t2d(const string&,const string&,const string&,int,int,int,double);
// get rupture from (displacement exceeding threshold)
void slip2rf(const string&,const string&,const string&,int,int,int,double,double,int,double);
// get total slip distribution at time t from posfiles & RAW_SERIES
void slip2d_total_rs(const string&, const string&, const  string&, const string&,  const string&, int,int,double,double,int);
// rupture from from velocity files
void vel2rf(const string&,const string&,const string&,int,int,int,int,double,double,int,double,double,int,double);
// moment rate function from displacement data  using RAW_SERIES and initial pos. file
void slip2momrate_rs(const string&, const string&, const  string&, const string&,  const string&, int,int,int,double,double,int);
#endif // __SLIP2VTK2D_H
