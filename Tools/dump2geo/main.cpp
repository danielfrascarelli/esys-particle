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

using std::cerr;
using std::endl;
using std::ostringstream;

// --- project includes ---
#include "frame_geo.h"

int main(int argc,char** argv)
{
  string infilename;
  string ofilename;
  string basefilename;
  string listfilename;
  string brklistfilename;
  bool with_rot=false;
  bool options_valid=true;
  bool with_list=false;
  bool with_brklist=false;
  bool remove_xbonds=false;
  double xbond_dist=0.0;
  int t0_frame=0;
  int tstart=0;
  int nt=0;
  int dt=0;
  int ret=0;

  // process args
  int args_read=1;  
  while(args_read<argc){
    string option=string(argv[args_read]);
    if(option=="-i"){
      if(argc>args_read){
        basefilename=string(argv[args_read+1]);
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-list"){
      if(argc>args_read){
        listfilename=string(argv[args_read+1]);
        with_list=true;
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-brklist"){
      if(argc>args_read){
        brklistfilename=string(argv[args_read+1]);
        with_brklist=true;
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
    }  else if(option=="-t0"){
      if(argc>args_read){
        t0_frame=atoi(argv[args_read+1]);
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-t"){
      if(argc>args_read+2){
        tstart=atoi(argv[args_read+1]);
        nt=atoi(argv[args_read+2]);
        dt=atoi(argv[args_read+3]);
        args_read+=4;
      } else {
        options_valid=false;
      }
    } else if(option=="-rxb"){
      if(argc>args_read+1){
        remove_xbonds=true;
        xbond_dist=atof(argv[args_read+1]);
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-rot") {
      with_rot=true;
      options_valid=true;
      args_read++;
    } else {
      cerr << "Unknown option " << option << endl;
      options_valid = false;
      break;
    }
  } 

  if(options_valid){
    for(int i=0;i<nt;i++){
      // make filename
      ostringstream filename;
      filename << basefilename << "_t=" << tstart+i*dt << "_0.txt";     // write frame
      if(with_rot){
        ostringstream brklistname;
        if(with_brklist){
          brklistname << brklistfilename << "_t=" << tstart+i*dt << "_0.txt";
        }
        do_single_frame_geo_r(filename.str(),ofilename,i+t0_frame,with_list,listfilename,remove_xbonds,xbond_dist,with_brklist,brklistname.str());
      } else {
        do_single_frame_geo(filename.str(),ofilename,i+t0_frame,listfilename);
      }
    }
    ret=0;
  } else {
    cerr << "Options invalid" << endl;
    ret=1;
  }

  return ret;
}
