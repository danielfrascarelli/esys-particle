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

#ifndef __FRAME_H
#define __FRAME_H

// -- project includes --
#include "colormap.h"

// -- STL includes --
#include <vector>
#include <string>
#include <map>

using std::string;
using std::vector;
using std::map;

void do_single_frame(const string&,const string&,int,const vector<string>&,int,const ColorMap*,int,Vec3,Vec3,bool,bool);
void do_single_frame_r(const string&,const string&,int,const vector<string>&,int,const ColorMap*,int,Vec3,Vec3,bool,bool,bool,map<int,Vec3>&);
map<int,Vec3> get_frame_disp_r(const string&);

#endif // __FRAME_H
