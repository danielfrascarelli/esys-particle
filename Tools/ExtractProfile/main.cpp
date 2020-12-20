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
#include <iostream>
#include <cstdlib>

using std::cerr;
using std::endl;
using std::string;

// --- project includes ---
#include "read.h"

int main(int argc,char** argv)
{
  string infilename;
  string infilename2;
  string outfilename;
  int ret=0;
  int mintag=0;
  bool options_valid=true;
  bool debug_on=false;
  bool grad_on=false;
  bool poros_on=false;
  bool diff_on=false;
  double xmin=0.0,xmax=0.0;
  double ymin=0.0,ymax=0.0;
  double zmin=0.0,zmax=0.0;
  double csize=0.0;
  bool dim3_on=false;
  int nbin=0;
  int dir=0;
  int udim=0, gdim=0;
 
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
    } else if(option=="-diff"){
      if(argc>args_read+1){
        infilename=string(argv[args_read+1]);
        infilename2=string(argv[args_read+2]);
        diff_on=true;
        args_read+=3;
      } else {
        options_valid=false;
      }
    } else if(option=="-o"){
      if(argc>args_read){
        outfilename=string(argv[args_read+1]);
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-debug"){
      if(argc>=args_read){
        debug_on=true;
      } else {
        options_valid=false;
      }
      args_read++;
    } else if(option=="-grad"){
      if(argc>=args_read){
        grad_on=true;
      } else {
        options_valid=false;
      }
      args_read++;
    }else if(option=="-b"){
      if(argc>args_read+2){
        ymin=atof(argv[args_read+1]);
        ymax=atof(argv[args_read+2]);
        nbin=atoi(argv[args_read+3]);
        args_read+=4;
      } else { 
        options_valid=false;
      }
    }else if(option=="-grad_t"){
      if(argc>args_read+1){
        udim=atoi(argv[args_read+1]);
        gdim=atoi(argv[args_read+2]);
        args_read+=3;
        grad_on=true;
      } else { 
        options_valid=false;
      }
    }else if(option=="-3d"){
      if(argc>args_read+6){
        xmin=atof(argv[args_read+1]);
        xmax=atof(argv[args_read+2]);
        ymin=atof(argv[args_read+3]);
        ymax=atof(argv[args_read+4]);
        zmin=atof(argv[args_read+5]);
        zmax=atof(argv[args_read+6]);
        csize=atof(argv[args_read+7]);
        dim3_on=true;
        args_read+=8;
      } else { 
        options_valid=false;
      }
    }else if(option=="-poros"){
      if(argc>args_read+6){
        xmin=atof(argv[args_read+1]);
        xmax=atof(argv[args_read+2]);
        ymin=atof(argv[args_read+3]);
        ymax=atof(argv[args_read+4]);
        zmin=atof(argv[args_read+5]);
        zmax=atof(argv[args_read+6]);
        csize=atof(argv[args_read+7]);
        poros_on=true;
        args_read+=8;
      } else { 
        options_valid=false;
      }
    }else if(option=="-dir"){
      if(argc>args_read){
        dir=atoi(argv[args_read+1]);
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-mintag"){
      if(argc>args_read){
        mintag=atoi(argv[args_read+1]);
        args_read+=2;
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
    if(dim3_on){
      read_and_write_disp_grid(infilename,outfilename,xmin,xmax,ymin,ymax,zmin,zmax,csize,grad_on,udim,gdim);
      ret=0;
    } else if (poros_on) {
      read_and_write_poros_grid(infilename,outfilename,xmin,xmax,ymin,ymax,zmin,zmax,csize);
      ret=0;
    } else if (diff_on) {
      read_and_write_profile_rel(infilename,infilename2,outfilename,ymin,ymax,nbin,debug_on,grad_on,dir);
    } else {
      read_and_write_profile_r(infilename,outfilename,ymin,ymax,nbin,debug_on,grad_on,dir,mintag);
      ret=0;
    }
  } else {
    cerr << "Options invalid" << endl;
    ret=1;
  }

  return ret;
}
