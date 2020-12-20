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

using std::string;
using std::cerr;
using std::endl;
using std::ostringstream;

// --- project includes ---
#include "slip2vtk2d.h"

/*
  use
  -rs POS_UP POS_DOWN:  use pos.files & RAW_SERIES files for disp
  -i DISP_UP DISP_DOWN: displacement files
  -o OUTPUT_FILE
  -ro OUTPUT_FILE raw output
  -shift SHIFT circular shift data if slip not centered (default 0)
  -vtk output vtk format
  -t T0 NT DT
  -x  X0 X1 NX
  -tl t-dimension of outfput grid (if -vtk)
  -r: output rate
  -td T0 T1 : total slip beween t0 and t1
  -rf get rupture fron (from displacements)
  -rfv THR OFS get rupture front from velocities
  -mr output moment rate file
  -ip STEPS interpolate steps between data points (only for raw output) default : 1 == no added points
*/
int main(int argc,char** argv)
{
  string ofilename;
  string basefilename_up;
  string basefilename_down;
  string posfilename_up;
  string posfilename_down;
  bool options_valid=true;
  bool rate=false;
  bool convert_to_vtk=false;
  bool single_point=false;
  bool tlset=false;
  bool get_rf=false;
  bool get_rf_vel=false;
  bool use_raw_series=false;
  bool output_raw=false;
  bool get_total_disp=false;
  bool get_momrate=false;
  int tstart=0;
  int nt=0;
  int dt=0;
  int tend=0;
  double x0=0.0, px0=0.0;
  double x1=0.0, px1=0.0;
  double tlen=0.0; // length of the resulting data in t-direction
  double ofs=0.0;
  double thr=0.0;
  int shift=0;
  int nx=0, npx=0;
  int ips=1;
  
  // process args
  int args_read=1;  
  while((args_read<argc) &&  options_valid ){
    string option=string(argv[args_read]);
    if(option=="-i"){
      if(argc>args_read){
        basefilename_up=string(argv[args_read+1]);
        basefilename_down=string(argv[args_read+2]);
        args_read+=3;
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
    } else if(option=="-ip"){
      if(argc>args_read){
        ips=atoi(argv[args_read+1]);
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-shift"){
      if(argc>args_read){
        shift=atoi(argv[args_read+1]);
        args_read+=2;
      } else {
        options_valid=false;
      }
    }else if(option=="-ro"){
      if(argc>args_read){
        ofilename=argv[args_read+1];
        output_raw=true;
        args_read+=2;
      } else {
        options_valid=false;
      }
    }else if(option=="-sp"){
      if(argc>args_read){
        x0=atof(argv[args_read+1]);
        single_point=true;
        args_read+=2;
      } else {
        options_valid=false;
      }
    }else if(option=="-rf"){
      if(argc>args_read){
        thr=atof(argv[args_read+1]);
        get_rf=true;
        args_read+=2;
      } else {
        options_valid=false;
      }
    }else if(option=="-rfv"){
      if(argc>args_read+1){
        thr=atof(argv[args_read+1]);
        ofs=atof(argv[args_read+2]);
        get_rf_vel=true;
        args_read+=3;
      } else {
        options_valid=false;
      }
    }else if(option=="-t"){
      if(argc>args_read+2){
        tstart=atoi(argv[args_read+1]);
        nt=atoi(argv[args_read+2]);
        dt=atoi(argv[args_read+3]);
        args_read+=4;
      } else {
        options_valid=false;
      } 
    } else if(option=="-tl"){
      if(argc>args_read){
        tlen=atof(argv[args_read+1]);
        tlset=true;
        args_read+=2;
      } else {
        options_valid=false;
      } 
    } else if(option=="-rs"){
      if(argc>args_read){
        posfilename_up=argv[args_read+1];
        posfilename_down=argv[args_read+2];
        use_raw_series=true;
        args_read+=3;
      } else {
        options_valid=false;
      } 
    } else if(option=="-td"){
      if(argc>args_read){
        tstart=atoi(argv[args_read+1]);
        tend=atoi(argv[args_read+2]);
        get_total_disp=true;
        args_read+=3;
      } else {
        options_valid=false;
      } 
    } else if(option=="-x"){
      if(argc>args_read+2){
        x0=atof(argv[args_read+1]);
        x1=atof(argv[args_read+2]);
        nx=atoi(argv[args_read+3]);
        args_read+=4;
      } else {
        options_valid=false;
      }
    } else if(option=="-px"){
      if(argc>args_read+2){
        px0=atof(argv[args_read+1]);
        px1=atof(argv[args_read+2]);
        npx=atoi(argv[args_read+3]);
        args_read+=4;
      } else {
        options_valid=false;
      }
    } else if (option=="-r"){
      rate=true;
      args_read++;
    } else if (option=="-mr"){
      get_momrate=true;
      args_read++;
    } else if (option=="-vtk"){
      convert_to_vtk=true;
      args_read++;
    } else {
      cerr << "Unknown option " << option << endl;
      options_valid = false;
      break;
    }
        
  }
  if(convert_to_vtk){
    if(!tlset) tlen=double(nt);
    if(rate){
      if(use_raw_series){
        slip2vtk2d_rate_rs(basefilename_up,basefilename_down,posfilename_up,posfilename_down,ofilename,tstart,nt,dt,tlen,x0,x1,nx);
      } else {
        slip2vtk2d_rate(basefilename_up,basefilename_down,ofilename,tstart,nt,dt,tlen,x0,x1,nx);      
      }
    } else {
      slip2vtk2d(basefilename_up,basefilename_down,ofilename,tstart,nt,dt,tlen,x0,x1,nx);
    }
  } else if (single_point){
    slip_x2slip_t2d(basefilename_up,basefilename_down,ofilename,tstart,nt,dt,x0);
  } else if (get_rf){
    slip2rf(basefilename_up,basefilename_down,ofilename,tstart,nt,dt,x0,x1,nx,thr);
  }else if (get_total_disp){
    slip2d_total_rs(basefilename_up,basefilename_down,posfilename_up,posfilename_down,ofilename,tstart,tend,x0,x1,nx);
  }else if (output_raw){
    if(use_raw_series){
      if(rate){
        slip2raw2d_rate_rs(basefilename_up,basefilename_down,posfilename_up,posfilename_down,ofilename,tstart,nt,dt,x0,x1,nx);
      } else   {
        slip2raw2d_rs(basefilename_up,basefilename_down,posfilename_up,posfilename_down,ofilename,tstart,nt,dt,x0,x1,nx,shift,ips);
      }
    } else {
      slip2raw2d(basefilename_up,basefilename_down,ofilename,tstart,nt,dt,x0,x1,nx);
    }
  } else if (get_rf_vel) {
    vel2rf(basefilename_up,basefilename_down,ofilename,tstart,nt,dt,ofs,x0,x1,nx,px0,px1,npx,thr);
  } else if (get_momrate) {
    slip2momrate_rs(basefilename_up,basefilename_down,posfilename_up,posfilename_down,ofilename,tstart,nt,dt,x0,x1,nx);
  } else {
    std::cout << "What do you want ??" << std::endl;
  } 
}
