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

#ifndef __FRACFRAME_H
#define __FRACFRAME_H

//--- system includes ---
#include <string>
#include <vector>
#include <set>
#include <utility>
#include <map>

using std::string;
using std::vector;
using std::set;
using std::pair;
using std::map;

//--- project includes ---
#include "Foundation/vec3.h"



class FracFrame
{
 private:
  // bond data

  struct bdata
  {
    int id1,id2,tag;
  };
  
  class cmp
  {
  public:
    bool operator ()(const bdata&, const bdata&);
  };

  set<bdata,FracFrame::cmp> m_bonds;
  map<int,Vec3> posmap;
  map<int,float> radmap;
  map<int,int> tagmap;

  // private helper functions
  int get_version(const string&);
  vector<string> get_filenames(const string&, int);

 public:
  
  // types
  struct fdata
  {
    Vec3 pos;
    Vec3 normal;
    double size;
    double dist;
    int id1, id2;
    int ptag1, ptag2;
    int tag;
  };

 
  FracFrame();

  void readFile(const string&);
  void readFileRot(const string&,bool);
  void readFileTagged(const string&,int);
  void readFileRotTagged(const string&,int);
  vector<fdata> getFrac(FracFrame&);
};

#endif // __FRACFRAME_H
