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

using std::cerr;
using std::endl;
using std::string;

// --- project includes ---
#include "DataExtractor.h"

int main(int argc,char** argv)
{
  string infilename;
  string outfilename;
  string dataname;
  int ret;
  int xdim=0, ydim=0, zdim=0;
  double gsize=0.0, srange=0.0;
  bool options_valid=true;
  bool debug_on=false;

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
        outfilename=string(argv[args_read+1]);
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-strain"){
      if(argc>args_read){
        srange=atof(argv[args_read+1]);
        dataname="strain";
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-grid"){
      if(argc>args_read+3){
        xdim=atoi(argv[args_read+1]);
        ydim=atoi(argv[args_read+2]);
        zdim=atoi(argv[args_read+3]);
        gsize=atof(argv[args_read+4]);    
        args_read+=5;
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
    } else {
      cerr << "Unknown option " << option << endl;
      options_valid = false;
      break;
    }
  } 

  if(options_valid){
    DataExtractor D(xdim,ydim,zdim,gsize,Vec3(0.0,0.0,0.0),Vec3(0.0,0.0,0.0));//modified (fluid contents)
    D.read(infilename);
    D.StrainToTensorData(srange);
    D.MaxShearToScalarData();
    D.writeScalarDataVtk(outfilename,dataname);
    ret=0;
  } else {
    cerr << "Options invalid" << endl;
    ret=1;
  }

  return ret;
}
