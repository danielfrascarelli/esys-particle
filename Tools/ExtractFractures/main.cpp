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

/*
  -i inputfile_prefix
  -o outputfile_name
  -ip output fractures at _initial_ particle positions


*/
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

// --- Project includes ---
#include "fracframe.h"
#include "fracwriter.h"
#include "FitPlane.h"
#include "frac_dist.h"
#include "Geometry/Plane3D.h"

int main(int argc,char** argv)
{
  bool options_valid=true;
  bool is_rot=false;
  bool is_tagged=false;
  bool do_fitplane=false;
  bool write_distrib=false;
  bool write_profile=false;
  bool write_list=false;
  bool write_format_vtk=false;
  bool write_ip=false;
  string ofilename;
  string basefilename;
  string listbasename;
  string prof_file;
  int t0_frame=0;
  int tstart=0;
  int nt=0;
  int dt=0;
  int tag=0;
  double dist_min=0.0;
  double dist_max=0.0;
  double ymin=0.0, ymax=0.0;
  int prof_nbin=0;
  int dist_nbin=0;
  string dist_filename;

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
        is_tagged=true;
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-t0"){
      if(argc>args_read){
        t0_frame=atoi(argv[args_read+1]);
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-dist"){
      if(argc>args_read+3){
        write_distrib=true;
        dist_filename=string(argv[args_read+1]);
        dist_min=atof(argv[args_read+2]);
        dist_max=atof(argv[args_read+3]);
        dist_nbin=atoi(argv[args_read+4]);
        args_read+=5;
      } else {
        options_valid=false;
      }
    } else if(option=="-prof"){
      if(argc>args_read+3){
        write_profile=true;
        prof_file=string(argv[args_read+1]);
        ymin=atof(argv[args_read+2]);
        ymax=atof(argv[args_read+3]);
        prof_nbin=atoi(argv[args_read+4]);
        args_read+=5;
      } else {
        options_valid=false;
      }
    } else if(option=="-list"){
      if(argc>args_read){
        listbasename=string(argv[args_read+1]);
        write_list=true;
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-rot"){
      is_rot=true;
      args_read++;
    } else if(option=="-pl"){
      do_fitplane=true;
      args_read++;
    } else if (option=="-vtk"){
      write_format_vtk=true;
      args_read++;
    } else if (option=="-ip"){
      write_ip=true;
      args_read++;
    } else if(option=="-t"){
      if(argc>args_read+2){
        tstart=atoi(argv[args_read+1]);
        nt=atoi(argv[args_read+2]);
        dt=atoi(argv[args_read+3]);
        args_read+=4;
        if(nt<2){
          options_valid=false;
        }
      } else {
        options_valid=false;
      }
    } else {
      cerr << "Unknown option " << option << endl;
      options_valid = false;
      break;
    }
  }

  // DO STUFF
  if(options_valid){
    FracWriter FW;
    FracFrame F1;
    ostringstream infilename1;

    // make filename
    infilename1 << basefilename << "_t=" << tstart << "_0.txt";
     
    if(is_tagged){
      if(is_rot){
        F1.readFileRotTagged(infilename1.str(),tag);
      } else {
        F1.readFileTagged(infilename1.str(),tag);
      }
    } else {
      if(is_rot){
        F1.readFileRot(infilename1.str(),write_ip);
      } else {
        F1.readFile(infilename1.str());
      }
    }
    for(int i=1;i<nt;i++){
      FracFrame F2;
      // make filename
      ostringstream infilename2;
      infilename2 << basefilename << "_t=" << tstart+i*dt << "_0.txt";
      if(is_tagged){
        if(is_rot){
          F2.readFileRotTagged(infilename2.str(),tag);
        } else {
          F2.readFileTagged(infilename2.str(),tag);
        }
      } else {
        if(is_rot){
          F2.readFileRot(infilename2.str(),write_ip);
        } else {
          F2.readFile(infilename2.str());
        }  
      }
      vector<FracFrame::fdata> newfrac=F1.getFrac(F2);
      F1=F2;
      if(write_distrib){
        FracDist fd=FracDist(newfrac,dist_min,dist_max,dist_nbin);
        fd.write(dist_filename);
      } else {
        FW.addData(newfrac,tstart+i*dt);
        if (write_list) {
          ostringstream listfilename;
          listfilename << listbasename << "_t=" << tstart+i*dt << "_0.txt";
          FW.writeParticleList(listfilename.str());
        } 
        if(do_fitplane){
          Plane3D Pl=fitPlaneToFracture(newfrac);
          cout << "Normal: [" << Pl.getNormal() << "] Position: [" << Pl.getPos() << "]" << endl; 
          cout << "Spanning vectors: [" << Pl.GetU() << "] [" << Pl.GetV() << "]" << endl; 
          FW.addPlane(Pl);
        }
      }
    }
    if(!write_distrib){
      if(write_profile){
        FW.writeProfile(ymin,ymax,prof_nbin,prof_file);
      } else if (write_format_vtk) {
        FW.write(ofilename);
      } else {
	FW.writeText(ofilename);
      }
    }
  }
  return 0;
}
