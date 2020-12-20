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
#include "gvmap.h"
#include "rextract.h"

int main(int argc,char** argv)
{
  string infilename;
  string ofilename;
  string basefilename;
//int t0_frame=0;
  int tstart=0;
  int nt=0;
  int dt=0;
  int ret=0;
  int tag=0;
  double x0=0.0,dx=0.0;
  int nb=0;
  bool options_valid=true;
  bool rex=false;
  bool rp=false;
  bool binned=false;

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
    } else if(option=="-o"){
      if(argc>args_read){
        ofilename=argv[args_read+1];
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-tag"){
      if(argc>args_read){
        tag=atoi(argv[args_read+1]);
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-r"){
      rex=true;
      args_read++;
    } else if(option=="-rl"){
      rp=true;
      args_read++;
    } else if(option=="-t"){
      if(argc>args_read+2){
        tstart=atoi(argv[args_read+1]);
        nt=atoi(argv[args_read+2]);
        dt=atoi(argv[args_read+3]);
        args_read+=4;
      } else {
        options_valid=false;
      }
    } else if(option=="-b"){
      if(argc>args_read+2){
        x0=atof(argv[args_read+1]);
        dx=atof(argv[args_read+2]);
        nb=atoi(argv[args_read+3]);
        binned=true;
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
    if(rp){ // extract from per-particle rotations
      Rextract R(basefilename,ofilename,tag);
      for(int i=0;i<nt;i++){
        R.read_frame(tstart+i*dt);
        if((i%10)==0) cout << "frame " << i << " [ " << tstart+i*dt << " ]" << endl;
      }
      if(binned){
        R.write_data_bin(x0,dx,nb);
      } else {
        R.write_data();
      }
    } else {
      GVMap gvm(tag);
      for(int i=0;i<nt;i++){
        // make filename
        ostringstream filename;
        filename << basefilename << "_t=" << tstart+i*dt << "_0.txt";
        gvm.read(filename.str(),rex);
        gvm.calc();
      }
    }
    ret=0;
  } else {
    ret=1;
  }

  return ret;
}
