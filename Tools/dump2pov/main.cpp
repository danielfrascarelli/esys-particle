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
#include <map>
#include <cstdlib>
#include <utility>
#include <vector>

// --- Project includes ---
#include "frame.h"
#include "mesh.h"
#include "colormap.h"
#include "colormap3.h"
#include "geocolormap.h"
#include "campos.h"
#include "../../Foundation/vec3.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ostringstream;
using std::atoi;
using std::vector;

int main(int argc,char** argv)
{
  vector<string> includefilename;
  string infilename;
  string ofilename;
  string basefilename;
  int t0_frame=0;
  int disp0_frame=0;
  int tstart=0;
  int nt=0;
  int dt=0;
  int camdir=0;
  bool options_valid=true;
  Vec3 campos,lookat;
  bool cam_set=false;
  bool with_mesh=false;
  bool with_rot=false;
  ColorMap* color_map=NULL;
  int field=-1;
  string camerafilename;
  bool with_campos=false;
  CameraPos* camera_descr=NULL;
  bool bonds_only=false;
  bool no_bonds=false;

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
    } else if(option=="-t0"){
      if(argc>args_read){
        t0_frame=atoi(argv[args_read+1]);
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-disp0"){
      if(argc>args_read){
        disp0_frame=atoi(argv[args_read+1]);
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-inc"){
      if(argc>args_read){
        includefilename.push_back(string(argv[args_read+1]));
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-campos"){
      if(argc>args_read){
        camerafilename=string(argv[args_read+1]);
        with_campos=true;
        cam_set=true;
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
    } else if(option=="-dir"){
      if(argc>args_read){
        camdir=atoi(argv[args_read+1]);
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-field"){
      if(argc>args_read){
        field=atoi(argv[args_read+1]);
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-camera"){
      if(argc>args_read+5){
        double lookx=atof(argv[args_read+1]);
        double looky=atof(argv[args_read+2]);
        double lookz=atof(argv[args_read+3]);
        double cposx=atof(argv[args_read+4]);
        double cposy=atof(argv[args_read+5]);
        double cposz=atof(argv[args_read+6]);
        cam_set=true;
        lookat=Vec3(lookx,looky,lookz);
        campos=Vec3(cposx,cposy,cposz);
        args_read+=7;
      } else {
        options_valid=false;
      }
    } else if(option=="-color"){
      if(argc>args_read+7){
        double r0=atof(argv[args_read+1]);
        double g0=atof(argv[args_read+2]);
        double b0=atof(argv[args_read+3]);
        double r1=atof(argv[args_read+4]);
        double g1=atof(argv[args_read+5]);
        double b1=atof(argv[args_read+6]);
        double x0=atof(argv[args_read+7]);
        double x1=atof(argv[args_read+8]);
        color_map=new ColorMap(Vec3(r0,g0,b0),Vec3(r1,g1,b1),x0,x1);
        args_read+=9;
      } else {
        options_valid=false;
      } 
    } else if(option=="-color3"){
      if(argc>args_read+11){
        double r0=atof(argv[args_read+1]);
        double g0=atof(argv[args_read+2]);
        double b0=atof(argv[args_read+3]);
        double r1=atof(argv[args_read+4]);
        double g1=atof(argv[args_read+5]);
        double b1=atof(argv[args_read+6]);
        double r2=atof(argv[args_read+7]);
        double g2=atof(argv[args_read+8]);
        double b2=atof(argv[args_read+9]);
        double x0=atof(argv[args_read+10]);
        double x1=atof(argv[args_read+11]);
        double x2=atof(argv[args_read+12]);
        color_map=new ColorMap3(Vec3(r0,g0,b0),Vec3(r1,g1,b1),Vec3(r2,g2,b2),x0,x1,x2);
        args_read+=13;
      } else {
        options_valid=false;
      } 
    } else if(option=="-geo"){
      if(argc>args_read+3){
        double x0=atof(argv[args_read+1]);
        double x1=atof(argv[args_read+2]);
        int nl=atoi(argv[args_read+3]);
        double rd=atof(argv[args_read+4]);
        color_map=new GeoColorMap(Vec3(0.64,0.48,0.395),Vec3(0.40,0.34,0.265),x0,x1,nl,rd);
        field=9;
        args_read+=5;
      } else {
        options_valid=false;
      }
    } else if(option=="-mesh") {
      with_mesh=true;
      args_read++;
    } else if(option=="-rot") {
      with_rot=true;
      args_read++;
    } else if(option=="-bonds") {
      bonds_only=true;
      args_read++;
    } else if(option=="-nobonds") {
      no_bonds=true;
      args_read++;
    } else {
      cerr << "Unknown option " << option << endl;
      options_valid = false;
      break;
    }
  }

  if(options_valid){
    map<int,Vec3> oldposmap;
    // read camera description file
    if(with_campos){
      camera_descr=new CameraPos(camerafilename);
    } 
    if(with_rot){ 
      ostringstream filename;
      filename << basefilename << "_t=" << disp0_frame << "_0.txt"; 
      oldposmap=get_frame_disp_r(filename.str());
    }
    for(int i=0;i<nt;i++){
      // make filename
      ostringstream filename;
      filename << basefilename << "_t=" << tstart+i*dt << "_0.txt"; 
      if(with_mesh){
        // setup filenames 
        // - evil hack - use ???_1.txt as source file -> change if mesh is distributed
        ostringstream mesh_infilename;
        mesh_infilename << basefilename << "_t=" << tstart+i*dt << "_1.txt"; 
        ostringstream mesh_outfilename;
        mesh_outfilename << "mesh_include_t=" << i+t0_frame << "_1.inc"; 
        // write mesh include
        do_mesh(mesh_infilename.str(),mesh_outfilename.str());
        // add mesh to include files
        includefilename.push_back(mesh_outfilename.str());
      }
      // setup camera position
      if(with_campos){
        pair<Vec3,Vec3> cp=camera_descr->getCamPos(tstart+i*dt);
        campos=cp.first;
        lookat=cp.second;
      }
      // write frame
      if(with_rot){
        do_single_frame_r(filename.str(),ofilename,i+t0_frame,includefilename,camdir,color_map,field,lookat,campos,cam_set,bonds_only,no_bonds,oldposmap);
      } else {
        do_single_frame(filename.str(),ofilename,i+t0_frame,includefilename,camdir,color_map,field,lookat,campos,cam_set,bonds_only);
      }
    }
  }

  if(color_map!=NULL) delete color_map;
  if(camera_descr!=NULL) delete camera_descr;
  return 0;
}

