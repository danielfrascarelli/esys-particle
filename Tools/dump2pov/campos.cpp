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

#include "campos.h"

// -- IO includes --
#include <fstream>
#include <iostream>

using std::ifstream;
using std::make_pair;

CameraPos::CameraPos(const string& filename)
{
  ifstream infile(filename.c_str());
  int npos;

  infile >> npos;
  for(int i=0;i<npos;i++){
    double x,y,z;
    int t;

    infile >> t;
    infile >> x >> y >> z;
    Vec3 pos=Vec3(x,y,z);
    infile >> x >> y >> z;
    Vec3 lookat=Vec3(x,y,z);
    m_posmap.insert(make_pair(t,make_pair(pos,lookat)));
  }
}

pair<Vec3,Vec3> CameraPos::getCamPos(int ts)
{
  Vec3 pos;
  Vec3 lookat;

  map<int,pair<Vec3,Vec3> >::iterator iter=m_posmap.upper_bound(ts);
  if(iter==m_posmap.begin()){
    pos=(iter->second).first;
    lookat=(iter->second).second;
  } else if(iter!=m_posmap.end()) {
    map<int,pair<Vec3,Vec3> >::iterator iter2=iter;
    iter--;
    Vec3 pos1=(iter->second).first;
    Vec3 lookat1=(iter->second).second;
    Vec3 pos2=(iter2->second).first;
    Vec3 lookat2=(iter2->second).second;
    double ts1=double(iter->first);
    double ts2=double(iter2->first);
    double frac=(double(ts)-ts1)/(ts2-ts1);
    pos=pos1+frac*(pos2-pos1);
    lookat=lookat1+frac*(lookat2-lookat1);
  } else {
    iter--; // blow up if map is empty 
    pos=(iter->second).first;
    lookat=(iter->second).second;
  }
  return make_pair(pos,lookat);
}
