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

using std::cout;
using std::cerr;
using std::endl;
using std::string;

// --- Project includes ---
#include "vvf.h"
#include "vtk.h"

int main(int argc,char** argv)
{
  bool options_valid=true;
  string ofilename;
  string ifilename;
  bool write_to_vtk=false;
  bool write_to_vvf=false;

  // process args
  int args_read=1;  
  while(args_read<argc){
    string option=string(argv[args_read]);
    if(option=="-i"){
      if(argc>args_read){
        ifilename=string(argv[args_read+1]);
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
    } else if(option=="-vvf"){
      write_to_vvf=true;
      args_read++;
    } else if(option=="-vtk"){
      write_to_vtk=true;
      args_read++;
    } else {
      cerr << "Unknown option " << option << endl;
      options_valid = false;
      break;
    }
  }


  // DO STUFF
  if(options_valid){
    if(write_to_vvf){
      convert_to_vvf(ifilename,ofilename);
    } else if(write_to_vtk){
      convert_to_vtk(ifilename,ofilename);
    }
  }  

  return 0;
}
