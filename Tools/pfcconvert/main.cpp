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

// --- System includes ---
#include <string>
#include <sstream>
#include <iostream>
#include <cstdlib>

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ostringstream;

// --- project includes ---
#include "pfcc.h"
#include "../../Foundation/vec3.h"

int main(int argc,char** argv)
{
  string infilename;
  string ofilename;
  string bbxfilename;
  int ret=0;
  Vec3 bbx_min,bbx_max;
  int cx=0,cy=0,cz=0;
  double scale=0.0;

  bool options_valid=true;


  // process args
  int args_read=1;  
  while(args_read<argc){
    string option=string(argv[args_read]);
    if(option=="-i"){
      if(argc>args_read){
        infilename=string(argv[args_read+1]);
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-o"){
      if(argc>args_read){
        ofilename=argv[args_read+1];
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-bbx"){
      if(argc>args_read){
        bbxfilename=argv[args_read+1];
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-scale"){
      if(argc>args_read){
        scale=atof(argv[args_read+1]);
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-circ"){
      if(argc>args_read+3){
        cx=atoi(argv[args_read+1]);
        cy=atoi(argv[args_read+2]);
        cz=atoi(argv[args_read+3]);
        args_read+=4;
      } else {
        options_valid=false;
      }
    } else {
      cerr << "Unknown option " << option << endl;
      options_valid = false;
      break;
    }
  }

  if(options_valid){
    pair<Vec3,Vec3> bbx=read_bbx(bbxfilename);
    cout << "bounding box: " << bbx.first << " , " << bbx.second << endl;
    pfc_convert(infilename,ofilename,bbx.first,bbx.second,cx,cy,cz,scale);
    ret=0;
  } else {
    ret=1;
  }

  return ret;
}
