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

#include "frame.h"

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <map>
#include <sstream>
#include <cstdlib>
#include <utility>
#include "vec3.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
using std::copy;
using std::istream_iterator;
using std::back_inserter;
using std::map;
using std::ostringstream;
using std::atoi;
using std::pair;
using std::make_pair;

int get_version(const string& infilename)
{
  string dummystring;
  int version;
  ifstream headerfile(infilename.c_str()); 
  // read token  
  headerfile >> dummystring;

  if(dummystring=="V"){ // if V -> new version 
    headerfile >> version ;
    cout << "version : " << version << endl;
  } else {
    cout << "pre- V.1 version" << endl;
    version=0;
  }
  headerfile.close();

  return version;
}


void do_single_frame(const string& infilename,const string& outfilename,int iframe, const vector<string>& incfile,int camdir,const ColorMap *cm, int field,Vec3 clook,Vec3 cpos,bool use_cam,bool bonds_only)
{
  map<int,Vec3> posmap;
  map<int,float> radmap;
  map<int,Vec3> colormap;

  vector<pair<int,int> > bonds;

  int version=get_version(infilename.c_str());
  ifstream headerfile(infilename.c_str()); 
  float dummy,xmax,ymax,zmax,xmin,ymin,zmin;
  string dummystring;
  vector<string> filenames;

  if(version==0){
    headerfile >> dummy >> dummy >> dummy;
    headerfile >> dummystring >> dummy;
  } else if ((version==1) || (version==2) || (version==3)) {
    headerfile >> dummystring >> dummy;
    headerfile >> dummy >> dummy >> dummy;
    headerfile >> dummystring >> dummy;
  } else {
    cerr << "unknown checkpoint version " << version << endl;
  }

  // get bounding box
  headerfile >> dummystring;
  headerfile >> xmin >> ymin >> zmin >> xmax >> ymax >> zmax ;

  // ignore periodic bdry
  headerfile >> dummystring >> dummy >> dummy >> dummy;

  // ignore dimension
  headerfile >> dummystring >> dummystring;

  // get file names
  copy(istream_iterator<string>(headerfile),istream_iterator<string>(),back_inserter(filenames));
  
  headerfile.close();
  // test 
  for(vector<string>::iterator iter=filenames.begin();
      iter!=filenames.end();
      iter++){
    ifstream datafile(iter->c_str());
    std::cerr << iter->c_str() << std::endl; 
    vector<float> pdata;

    // get particles
    int npart;
    Vec3 pos;
    Vec3 oldpos;
    Vec3 initpos;
    Vec3 vel;
    Vec3 force;
    Vec3 circshift;
    double rad;
    double mass;
    int id;
    int tag;
    datafile >> npart;
    for(int i=0;i<npart;i++){
      datafile >> pos  >> rad >> id >> tag >> mass >> initpos >> oldpos >> vel >> force;
      if(version>1) datafile >> circshift ;
      posmap[id]=pos;
      radmap[id]=rad;
      double c_temp;
      switch(field){
      case 1: c_temp=pos.X(); break;
      case 2: c_temp=pos.Y(); break;
      case 3: c_temp=pos.Z(); break;
      case 4: c_temp=rad; break;
      case 5: c_temp=double(id); break;
      case 6: c_temp=double(tag); break;
      case 7: c_temp=mass; break;
      case 8: c_temp=initpos.X(); break;
      case 9: c_temp=initpos.Y(); break;
      case 10: c_temp=initpos.Z(); break;
      case 11: c_temp=vel.X(); break;
      case 12: c_temp=vel.Y(); break;
      case 13: c_temp=vel.Z(); break;
      case 14: c_temp=(pos-initpos).norm();; break;
      default: c_temp=0.0;
      }
      colormap[id]=cm->getColor(c_temp);
    }
      
    // get bonds
    int ngrp;
    int nbond;
    string type;

    datafile >> ngrp;
    if(version>0){
      datafile >> type;
      if(type=="Bonded"){
	datafile >> nbond;
	for(int i=0;i<nbond;i++){
	  int b1,b2;
	  
	  datafile >> b1 >> b2 >> dummy;
	  bonds.push_back(make_pair(b1,b2));
	}
      }
    } else { // pre - V1 snapshot -> assume bondend pair IG
      datafile >> nbond;
      for(int i=0;i<nbond;i++){
	int b1,b2;
	
	datafile >> b1 >> b2 >> dummy;
	bonds.push_back(make_pair(b1,b2));
      }
    }



    datafile.close();
  }
  
  // calculate camera position
  Vec3 lookat; 
  Vec3 loc;
  Vec3 light1_loc,light2_loc;

  if(use_cam){
    lookat=clook;
    loc=cpos;
    light1_loc=clook+Vec3(0.0,50.0,0.0);
    light1_loc=cpos+Vec3(0.0,10.0,0.0);
  } else {
    float cdist=0.0;
    switch(camdir){
    case 1: {// camera x
      loc=lookat+Vec3(xmax+cdist,0.0,0.0); 
      light1_loc=lookat+Vec3(xmax+2*cdist,0.0,cdist);
      light2_loc=lookat+Vec3(xmax+2*cdist,0.0,-1.0*cdist);
    } break; 
    case 2: {// camera y
      loc=lookat+Vec3(0.0,ymax+cdist,0.0); 
      light1_loc=lookat+Vec3(0.0,ymax+2*cdist,cdist);
      light2_loc=lookat+Vec3(0.0,ymax+2*cdist,-1.0*cdist);
    } break; 
    default: {// camera z
      loc=lookat+Vec3(0.0,0.0,zmax+cdist); 
      light1_loc=lookat+Vec3(cdist,0.0,zmax+2*cdist);
      light2_loc=lookat+Vec3(-1.0*cdist,0.0,zmax+2*cdist);
    } break; 
    }
  }
  // open output file
  ostringstream povfilename;
  povfilename << outfilename << iframe << ".pov";
  ofstream povfile(povfilename.str().c_str());
  // comment
  povfile << "//generated by dump2pov\n\n";
  // write camera, lights
  povfile << "#include \"colors.inc\" \n";
  for(vector<string>::const_iterator iter=incfile.begin();
      iter!=incfile.end();
      iter++){
    povfile << "#include \"" << *iter <<"\" \n";
  }
  povfile << "background { color White } \n";
  povfile << "camera { \n location <" << loc.X() << "," << loc.Y() << "," << loc.Z() << "> \n look_at  <" << lookat.X() << "," << lookat.Y() << "," << lookat.Z() << "> \n}\n";
  povfile << "light_source { <" << light1_loc.X() << "," << light1_loc.Y() << "," << light1_loc.Z() << "> color White}\n";
  povfile << "light_source { <" << light2_loc.X() << "," << light2_loc.Y() << "," << light2_loc.Z() << "> color White}\n\n";
  
  // write particles
  povfile << "union{\n";
  for(map<int,Vec3>::iterator iter=posmap.begin();
      iter!=posmap.end();
      iter++){
    Vec3 pos=iter->second;
    float rad=radmap[iter->first];
    Vec3 color=colormap[iter->first];
    povfile << "sphere { <" << pos.X() << " , " << pos.Y() << " , " << pos.Z() << ">, " << rad << "  " ;
    povfile << "pigment { rgb < " << color.X() << " , " << color.Y() << " , " << color.Z() << "> } }\n"; 
  }
  // write bonds
  for(vector<pair<int,int> >::iterator iter=bonds.begin();
      iter!=bonds.end();
      iter++){
    // get endpoints
    Vec3 p1=posmap[iter->first];
    Vec3 p2=posmap[iter->second];
    // get radii
    float r1=radmap[iter->first];
    float r2=radmap[iter->second];
    // sanity check
    if((p1-p2).norm()<2.0*(r1+r2)){ // cull circular bonds
      // calc color
      Vec3 color=0.5*(colormap[iter->first]+colormap[iter->second]);
      // write
      povfile << "cone { < " << p1.X() << "," << p1.Y() << "," << p1.Z() << "> " << r1*0.75 << " < " << p2.X() << "," << p2.Y() << "," << p2.Z() << "> " << r2*0.75 << " ";
      povfile << "pigment { rgb < " << color.X() << " , " << color.Y() << " , " << color.Z() << "> } }\n"; 
    }
  }
  povfile << "}\n";
  povfile.close();
}

void do_single_frame_r(const string& infilename,const string& outfilename,int iframe, const vector<string>& incfile,int camdir,const ColorMap *cm, int field,Vec3 clook,Vec3 cpos,bool use_cam,bool bonds_only,bool no_bonds,map<int,Vec3>& p_orig)
{
  map<int,Vec3> posmap;
  map<int,float> radmap;
  map<int,Vec3> colormap;

  vector<pair<int,int> > bonds;

  int version=get_version(infilename.c_str());
  ifstream headerfile(infilename.c_str()); 
  float dummy,xmax,ymax,zmax,xmin,ymin,zmin;
  string dummystring;
  vector<string> filenames;

  if(version==0){
    headerfile >> dummy >> dummy >> dummy;
    headerfile >> dummystring >> dummy; 
  } else if ((version==1) || (version==2) || (version==3)) {
    headerfile >> dummystring >> dummy;
    headerfile >> dummy >> dummy >> dummy;
    headerfile >> dummystring >> dummy;
  } else {
    cerr << "unknown checkpoint version " << version << endl;
  }

  // get bounding box
  headerfile >> dummystring;
  headerfile >> xmin >> ymin >> zmin >> xmax >> ymax >> zmax ;
  
  // ignore periodic bdry
  headerfile >> dummystring >> dummy >> dummy >> dummy;
  
  // ignore dimension
  headerfile >> dummystring >> dummystring;
  
  // get file names
  copy(istream_iterator<string>(headerfile),istream_iterator<string>(),back_inserter(filenames));
  
  headerfile.close();
  // test 
  double ct_min=1000;
  double ct_max=-1000;
  for(vector<string>::iterator iter=filenames.begin();
      iter!=filenames.end();
      iter++){
    ifstream datafile(iter->c_str());
    vector<float> pdata;
    
    // get particles
    int npart;
    Vec3 pos;
    Vec3 oldpos;
    Vec3 initpos;
    Vec3 vel;
    Vec3 force;
    Vec3 angvel;
    double rad;
    double mass;
    double q1,q2,q3,q4;
    int id;
    int tag;
    datafile >> npart;
    if(version<2){
      for(int i=0;i<npart;i++){
	datafile >> pos  >> rad >> id >> tag >> mass >> initpos >> oldpos >> vel >> force >> q1 >> q2 >> q3 >> q4 >> angvel;
	posmap[id]=pos;
	radmap[id]=rad;
	double c_temp;
	switch(field){
	case 1: c_temp=pos.X(); break;
	case 2: c_temp=pos.Y(); break;
	case 3: c_temp=pos.Z(); break;
	case 4: c_temp=rad; break;
	case 5: c_temp=double(id); break;
	case 6: c_temp=double(tag); break;
	case 7: c_temp=mass; break;
	case 8: c_temp=initpos.X(); break;
	case 9: c_temp=initpos.Y(); break;
	case 10: c_temp=initpos.Z(); break;
	case 11: c_temp=vel.X(); break;
	case 12: c_temp=vel.Y(); break;
	case 13: c_temp=vel.Z(); break;
	case 14: c_temp=(pos-initpos).norm(); break;
	case 15: c_temp=(pos-p_orig[id]).norm(); break;
	case 16: c_temp=p_orig[id].Z(); break;
	default: c_temp=0.0;
	}
	ct_max=c_temp > ct_max ? c_temp : ct_max;
	ct_min=c_temp < ct_min ? c_temp : ct_min;
	colormap[id]=cm->getColor(c_temp);
      }
    } else {
      Vec3 circ_shift;
      for(int i=0;i<npart;i++){
	datafile >> pos  >> rad >> id >> tag >> mass >> initpos >> oldpos >> vel >> force >> circ_shift >> q1 >> q2 >> q3 >> q4 >> angvel;
	posmap[id]=pos;
	radmap[id]=rad;
	double c_temp;
	switch(field){
	case 1: c_temp=pos.X(); break;
	case 2: c_temp=pos.Y(); break;
	case 3: c_temp=pos.Z(); break;
	case 4: c_temp=rad; break;
	case 5: c_temp=double(id); break;
	case 6: c_temp=double(tag); break;
	case 7: c_temp=mass; break;
	case 8: c_temp=initpos.X(); break;
	case 9: c_temp=initpos.Y(); break;
	case 10: c_temp=initpos.Z(); break;
	case 11: c_temp=vel.X(); break;
	case 12: c_temp=vel.Y(); break;
	case 13: c_temp=vel.Z(); break;
	case 14: c_temp=(pos-initpos).norm(); break;
	case 15: c_temp=(pos-p_orig[id]).norm(); break;
	case 16: c_temp=p_orig[id].Z(); break;
	default: c_temp=0.0;
	}
	ct_max=c_temp > ct_max ? c_temp : ct_max;
	ct_min=c_temp < ct_min ? c_temp : ct_min;
	colormap[id]=cm->getColor(c_temp);
      }
    }
    std::cout << "min/max : " << ct_min << " / " << ct_max << std::endl;

    // get bonds
    int ngrp;
    int nbond;
    string type;

    datafile >> ngrp;
    if(version>0){
      datafile >> type;
      if((type=="Bonded") || (type=="RotBonded")){
	datafile >> nbond;
	for(int i=0;i<nbond;i++){
	  int b1,b2;
	  
	  datafile >> b1 >> b2 >> dummy;
	  bonds.push_back(make_pair(b1,b2));
	}
      }
    } else { // pre - V1 snapshot -> assume bondend pair IG
      datafile >> nbond;
      for(int i=0;i<nbond;i++){
	int b1,b2;
	
	datafile >> b1 >> b2 >> dummy;
	bonds.push_back(make_pair(b1,b2));
      }
    }

    datafile.close();
  }
  
  // calculate camera position
  Vec3 lookat; 
  Vec3 loc;
  Vec3 light1_loc,light2_loc;

  if(use_cam){
    lookat=clook;
    loc=cpos;
    light1_loc=clook+Vec3(0.0,50.0,0.0);
    light1_loc=cpos+Vec3(0.0,10.0,0.0);
  } else {
    float cdist=0.0;
    switch(camdir){
    case 1: {// camera x
      loc=lookat+Vec3(xmax+cdist,0.0,0.0); 
      light1_loc=lookat+Vec3(xmax+2*cdist,0.0,cdist);
      light2_loc=lookat+Vec3(xmax+2*cdist,0.0,-1.0*cdist);
    } break; 
    case 2: {// camera y
      loc=lookat+Vec3(0.0,ymax+cdist,0.0); 
      light1_loc=lookat+Vec3(0.0,ymax+2*cdist,cdist);
      light2_loc=lookat+Vec3(0.0,ymax+2*cdist,-1.0*cdist);
    } break; 
    default: {// camera z
      loc=lookat+Vec3(0.0,0.0,zmax+cdist); 
      light1_loc=lookat+Vec3(cdist,0.0,zmax+2*cdist);
      light2_loc=lookat+Vec3(-1.0*cdist,0.0,zmax+2*cdist);
    } break; 
    }
  }
  // open output file
  ostringstream povfilename;
  povfilename << outfilename << iframe << ".pov";
  ofstream povfile(povfilename.str().c_str());
  // comment
  povfile << "//generated by dump2pov\n\n";
  // write camera, lights
  povfile << "#include \"colors.inc\" \n";
  for(vector<string>::const_iterator iter=incfile.begin();
      iter!=incfile.end();
      iter++){
    povfile << "#include \"" << *iter <<"\" \n";
  }
  povfile << "background { color White } \n";
 //  povfile << "camera { \n location <" << loc.X() << "," << loc.Y() << "," << loc.Z() << "> \n look_at  <" << lookat.X() << "," << lookat.Y() << "," << lookat.Z() << "> \n}\n";
//   povfile << "light_source { <" << light1_loc.X() << "," << light1_loc.Y() << "," << light1_loc.Z() << "> color White}\n";
//   povfile << "light_source { <" << light2_loc.X() << "," << light2_loc.Y() << "," << light2_loc.Z() << "> color White}\n\n";
  
  // write particles
  povfile << "union{\n";
  if(!bonds_only){
    for(map<int,Vec3>::iterator iter=posmap.begin();
	iter!=posmap.end();
	iter++){
      Vec3 pos=iter->second;
      float rad=radmap[iter->first];
      Vec3 color=colormap[iter->first];
      povfile << "sphere { <" << pos.X() << " , " << pos.Y() << " , " << pos.Z() << ">, " << rad << "  " ;
      povfile << "pigment { rgb < " << color.X() << " , " << color.Y() << " , " << color.Z() << "> } }\n"; 
    }
  }
  // write bonds
  if(!no_bonds){
    for(vector<pair<int,int> >::iterator iter=bonds.begin();
	iter!=bonds.end();
	iter++){
      // get endpoints
      Vec3 p1=posmap[iter->first];
      Vec3 p2=posmap[iter->second];
      // get radii
      float r1=radmap[iter->first];
      float r2=radmap[iter->second];
      // sanity check
      if((p1-p2).norm()<2.0*(r1+r2)){ // cull circular bonds
	// calc color
	Vec3 color=0.5*(colormap[iter->first]+colormap[iter->second]);
	// write
	povfile << "cone { < " << p1.X() << "," << p1.Y() << "," << p1.Z() << "> " << r1*0.75 << " < " << p2.X() << "," << p2.Y() << "," << p2.Z() << "> " << r2*0.75 << " ";
	povfile << "pigment { rgb < " << color.X() << " , " << color.Y() << " , " << color.Z() << "> } }\n"; 
      }
    }
  }
  povfile << "}\n";
  povfile.close();
}

map<int,Vec3> get_frame_disp_r(const string& infilename)
{
  map<int,Vec3> posmap;

  int version=get_version(infilename.c_str());
  ifstream headerfile(infilename.c_str()); 
  float dummy,xmax,ymax,zmax,xmin,ymin,zmin;
  string dummystring;
  vector<string> filenames;

  if(version==0){
    headerfile >> dummy >> dummy >> dummy;
    headerfile >> dummystring >> dummy; 
  } else if ((version==1) || (version==2)|| (version==3)) {
    headerfile >> dummystring >> dummy;
    headerfile >> dummy >> dummy >> dummy;
    headerfile >> dummystring >> dummy;
  } else {
    cerr << "unknown checkpoint version " << version << endl;
  }

  // get bounding box
  headerfile >> dummystring;
  headerfile >> xmin >> ymin >> zmin >> xmax >> ymax >> zmax ;
  
  // ignore periodic bdry
  headerfile >> dummystring >> dummy >> dummy >> dummy;
  
  // ignore dimension
  headerfile >> dummystring >> dummystring;
  
  // get file names
  copy(istream_iterator<string>(headerfile),istream_iterator<string>(),back_inserter(filenames));
  
  headerfile.close();
  for(vector<string>::iterator iter=filenames.begin();
      iter!=filenames.end();
      iter++){
    ifstream datafile(iter->c_str());
    vector<float> pdata;
    
    // get particles
    int npart;
    Vec3 pos;
    Vec3 oldpos;
    Vec3 initpos;
    Vec3 vel;
    Vec3 force;
    Vec3 angvel;
    double rad;
    double mass;
    double q1,q2,q3,q4;
    int id;
    int tag;
    datafile >> npart;
    if(version<2){
      for(int i=0;i<npart;i++){
	datafile >> pos  >> rad >> id >> tag >> mass >> initpos >> oldpos >> vel >> force >> q1 >> q2 >> q3 >> q4 >> angvel;
	posmap[id]=pos;
      }
    } else {
      Vec3 circ_shift;
      for(int i=0;i<npart;i++){
	datafile >> pos  >> rad >> id >> tag >> mass >> initpos >> oldpos >> vel >> force >> circ_shift >> q1 >> q2 >> q3 >> q4 >> angvel;
	posmap[id]=pos;
      }
    }
    datafile.close();
  }
  
  return posmap;
}
