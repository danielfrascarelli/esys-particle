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
frame_geo.cpp:
  Written by Feng Chen, Vince Boros and Michele Griffa between
  November 2010 and January 2011, based on ../dump2vtk/frame_vtk.cpp
  (https://answers.launchpad.net/esys-particle/+question/132497).
*/

#include "frame_geo.h"

// --- project includes ---
#include "vec3.h"


// --- STL includes ---
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <map>
#include <utility>
#include <sstream>

using std::vector;
using std::string;
using std::ifstream;
using std::istream_iterator;
using std::back_inserter;
using std::map;
using std::pair;
using std::make_pair;
using std::ostringstream;

string dimension;

struct nr_part{
  Vec3 pos;
  Vec3 vel;
  Vec3 force;
  double rad;
  double mass;
  int tag;
  int id;
};

struct r_part{
  Vec3 pos;
  Vec3 init_pos;
  Vec3 vel;
  Vec3 force;
  Vec3 circ_shift;
  double rad;
  double mass;
  int tag;
  int id;
  double q1,q2,q3,q4;
  Vec3 angvel;
  int proc_id;
};

struct bond
{
  int id1;
  int id2;
  int tag;

  bond(int,int,int);
};

bond::bond(int i1,int i2,int t)
{
  id1=i1;
  id2=i2;
  tag=t;
}

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

vector<string> get_filenames(const string& infilename, int version, string& dimension, float* geo_pbdry=NULL, double* bdbx=NULL)
{
  cout << "infilename : " << infilename << endl;
  ifstream headerfile(infilename.c_str()); 
  float dummy,xmax,ymax,zmax,xmin,ymin,zmin,geo_version;
  vector<string> filenames;
  string dummystring;

  if(version==0){
    headerfile >> dummy >> dummy >> dummy;
    headerfile >> dummystring >> dummy;
  } else if ((version==1)||(version==2) || (version==3)) {
    headerfile >> dummystring >> dummy;
    headerfile >> dummy >> dummy >> dummy;
    headerfile >> dummystring >> geo_version;
  } else {
    cerr << "unknown checkpoint version " << version << endl;
  }
  // get bounding box
  headerfile >> dummystring;
  headerfile >> xmin >> ymin >> zmin >> xmax >> ymax >> zmax ;
  
  //double boundingBox[6]={xmin,ymin,zmin,xmax,ymax,zmax};
  if (bdbx!=NULL)
  {
     bdbx[0]=xmin,bdbx[1]=ymin,bdbx[2]=zmin;
     bdbx[3]=xmax,bdbx[4]=ymax,bdbx[5]=zmax;
  }
  //std::cerr << bdbx[0] << " " << bdbx[1]<< " " << bdbx[2]<< " " << bdbx[3]<< " " << bdbx[4]<< " " << bdbx[5] << endl;

  //geo_pbdry[0] = geo_version;
  geo_pbdry[0] = 1.2;

  // periodic bdry
  headerfile >> dummystring >> geo_pbdry[1] >> geo_pbdry[2] >> geo_pbdry[3];

  // v. 1.1 geometry files didn't have dimension info
  if(geo_version>1.15){
    headerfile >> dummystring >> dimension;
  }

  // get file names
  copy(istream_iterator<string>(headerfile),istream_iterator<string>(),back_inserter(filenames));
  
  headerfile.close();
  
  //==================================
  
  //==================================

  cout << "nr. of filenames: " << filenames.size() << endl;
  return filenames;
}

void do_single_frame_geo(const string& infilename,const string& outfilename,int iframe,const string& listfilename)
{
  int version=get_version(infilename);
  string dimension;
  float geo_pbdry[4];
  double bdbx[6];
  vector<string> filenames=get_filenames(infilename,version,dimension,geo_pbdry,bdbx);
  
  std::cerr << bdbx[0] << " " << bdbx[1]<< " " << bdbx[2]<< " " << bdbx[3]<< " " << bdbx[4]<< " " << bdbx[5] << endl;
  
  map<int,nr_part> datamap;
  vector<bond> bonds;
  map<int,int> id2idx;
  map<int,int> idx2id;

  for(vector<string>::iterator iter=filenames.begin();
      iter!=filenames.end();
      iter++){
    cout << *iter << endl;
    ifstream datafile(iter->c_str());
    vector<double> pdata;
    // get particles
    Vec3 circ_shift;
    int npart;
    nr_part data;
    int id;
    Vec3 pos;
    Vec3 oldpos;
    Vec3 initpos;
    datafile >> npart;
    if(version<2){
      for(int i=0;i<npart;i++){
        datafile >> data.pos  >> data.rad >> id >> data.tag >> data.mass >> initpos >> oldpos >> data.vel >> data.force;
        data.id = id;
        datamap[id]=data;
      }
    } else {
      for(int i=0;i<npart;i++){
        datafile >> data.pos  >> data.rad >> id >> data.tag >> data.mass >> initpos >> oldpos >> data.vel >> data.force >> circ_shift;
        data.id = id;
        datamap[id]=data;
      }
    }

    // get bonds
    int ngrp;
    int nbond;
    string type;

    datafile >> ngrp;
    for(int ni=0;ni<ngrp;ni++){
      if(version>0){
        datafile >> type;
        if(type=="Bonded"){
          datafile >> nbond;
          for(int i=0;i<nbond;i++){
            int b1,b2,tag;
            datafile >> b1 >> b2 >> tag;
            bonds.push_back(bond(b1,b2,tag));
          }
        }
      } else { // pre - V1 snapshot -> assume bondend pair IG
        datafile >> nbond;
        for(int i=0;i<nbond;i++){
          int b1,b2,tag;
          datafile >> b1 >> b2 >> tag;
          bonds.push_back(bond(b1,b2,tag));
        }
      }
    }
    datafile.close();
  }
  // generate reverse mapping between particle id and point idx
  int count=0;
  for(map<int,nr_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    id2idx.insert(make_pair(iter->first,count));
    idx2id.insert(make_pair(count, iter->first));
    count++;
  }
  // open output file
  ostringstream geofilename;
  geofilename << outfilename << iframe << ".geo";
  ofstream geofile(geofilename.str().c_str());
  geofile << "LSMGeometry " << geo_pbdry[0] << endl;
  geofile << "BoundingBox " << bdbx[0] << " " << bdbx[1] << " " << bdbx[2] << " " << bdbx[3] << " " << bdbx[4] << " " << bdbx[5] << endl;
  geofile << "PeriodicBoundaries " << geo_pbdry[1] << " " << geo_pbdry[2] << " " << geo_pbdry[3] << endl;
  geofile << "Dimension " << dimension << endl;
  geofile << "BeginParticles" << endl;
  geofile << "Simple" << endl;
  geofile << datamap.size() << endl;
  for(map<int,nr_part>::iterator iter=datamap.begin(); iter!=datamap.end(); iter++){
    //geofile << (iter->second).pos << " "<< (iter->second).rad << " "<< geopid++ << " " << (iter->second).tag << endl;
    geofile << (iter->second).pos << " "<< (iter->second).rad << " "<< (iter->second).id << " " << (iter->second).tag << endl;
//    geopid++;
  }
  geofile << "EndParticles" << endl;
  //write bond info
  geofile << "BeginConnect" << endl;
  geofile << bonds.size() << endl;
  for(vector<bond>::iterator iter=bonds.begin(); iter!=bonds.end(); iter++){
    geofile << iter->id1 << " " << iter->id2 << " " << iter->tag << endl;
  }                                 
  geofile << "EndConnect" << endl;
  geofile.close();
}

void do_single_frame_geo_r(const string& infilename,const string& outfilename,int iframe,bool with_list,const string& listfilename,bool remove_xbonds,double bond_remove_dist,bool with_brklist,const string& brklistname)
{  
  int version=get_version(infilename);
  string dimension;
  float geo_pbdry[4];
  double bdbx[6];
  vector<string> filenames=get_filenames(infilename,version,dimension,geo_pbdry,bdbx);
  
  std::cerr << bdbx[0] << " " << bdbx[1]<< " " << bdbx[2]<< " " << bdbx[3]<< " " << bdbx[4]<< " " << bdbx[5] << endl;
  
  map<int,r_part> datamap;
  map<int,int> ng_id_map;
  map<int,float> ng_mass_map;
  vector<bond> bonds;
  map<int,int> id2idx;
  map<int,int> id_brk_map;

  int proc_cnt=0;
  for(vector<string>::iterator iter=filenames.begin();
      iter!=filenames.end();
      iter++){
    proc_cnt++;
    cout << *iter << " , " << proc_cnt << endl;
    ifstream datafile(iter->c_str());
    vector<double> pdata;
    // get particles
    int npart;
    r_part data;
    int id;
    Vec3 pos;
    Vec3 oldpos;
    datafile >> npart;
    if(version < 2){
      for(int i=0;i<npart;i++){
        datafile >> data.pos  >> data.rad >> data.id >> data.tag >> data.mass >> data.init_pos >> oldpos >> data.vel >> data.force
          >> data.q1 >> data.q2 >> data.q3 >> data.q4 >> data.angvel;
        data.proc_id=proc_cnt;
        id = data.id;
        datamap[id]=data;
      } 
    } else {
     for(int i=0;i<npart;i++){
       datafile >> data.pos  >> data.rad >> data.id >> data.tag >> data.mass >> data.init_pos >> oldpos >> data.vel >> data.force >> data.circ_shift
        >> data.q1 >> data.q2 >> data.q3 >> data.q4 >> data.angvel;
       data.proc_id=proc_cnt;
       id = data.id;
       datamap[id]=data;
      } 
    }  

    // get bonds
    int ngrp;
    int nbond;
    string type;

    
    datafile >> ngrp;
    for(int ni=0;ni<ngrp;ni++){
      if(version>0){
        datafile >> type;
        if((type=="Bonded") || (type=="RotBonded")){
          datafile >> nbond;
          if(remove_xbonds){
            for(int i=0;i<nbond;i++){
              int b1,b2,dummy;
              datafile >> b1 >> b2 >> dummy;
              double dist=(datamap[b1].pos-datamap[b2].pos).norm();
              if(dist < bond_remove_dist){
                bonds.push_back(bond(b1,b2,ni));
              }
            }
          } else {
            for(int i=0;i<nbond;i++){
              int b1,b2,dummy;
              datafile >> b1 >> b2 >> dummy;
              bonds.push_back(bond(b1,b2,ni));
            }
          }
        }
      } else { // pre - V1 snapshot -> assume bondend pair IG
        datafile >> nbond;
        for(int i=0;i<nbond;i++){
          int b1,b2,tag;
          datafile >> b1 >> b2 >> tag;
          bonds.push_back(bond(b1,b2,tag));
        }
      }
    }
    datafile.close();
  }

  // generate reverse mapping between particle id and point idx
  int count=0;
  for(map<int,r_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    id2idx.insert(make_pair(iter->first,count));
    count++;
  }

  // if option set, read list file into ng_id_map;
  if(with_list){
    ifstream listfile(listfilename.c_str());
    int id, nid;
    float gmass;
    while(!listfile.eof()){
      listfile >> id >> nid >> gmass;
      ng_id_map[id]=nid;
      ng_mass_map[id]=gmass;
    }
    listfile.close();
  }
  // if option set, read fracture list file int id_brk_map;
  if(with_brklist){
    std::cerr << "opening " << brklistname << std::endl;
    ifstream brklistfile(brklistname.c_str());
    int id, nbrk;
    while(!brklistfile.eof()){
      brklistfile >> id >> nbrk;
      id_brk_map[id]=nbrk;
    }
    brklistfile.close();
  }
  
  //std::cerr << "I am being called!!!!!!!!!!!!!!!!!!" << std::endl;
  
  // open output file
  ostringstream geofilename;
  geofilename << outfilename << iframe << ".geo";
  ofstream geofile(geofilename.str().c_str());
  geofile << "LSMGeometry " << geo_pbdry[0] << endl;
  geofile << "BoundingBox " << bdbx[0] << " " << bdbx[1] << " " << bdbx[2] << " " << bdbx[3] << " " << bdbx[4] << " " << bdbx[5] << endl;
  geofile << "PeriodicBoundaries " << geo_pbdry[1] << " " << geo_pbdry[2] << " " << geo_pbdry[3] << endl;
  geofile << "Dimension " << dimension << endl;
  geofile << "BeginParticles" << endl;
  geofile << "Simple" << endl;
  geofile << datamap.size() << endl;
  for(map<int,r_part>::iterator iter=datamap.begin(); iter!=datamap.end(); iter++){
    //geofile << (iter->second).pos << " "<< (iter->second).rad << " "<< geopid++ << " " << (iter->second).tag << endl;
    geofile << (iter->second).pos << " "<< (iter->second).rad << " "<< (iter->second).id << " " << (iter->second).tag << endl;
  }
  geofile << "EndParticles" << endl;
  //write bond info
  geofile << "BeginConnect" << endl;
  geofile << bonds.size() << endl;
  for(vector<bond>::iterator iter=bonds.begin(); iter!=bonds.end(); iter++){
    geofile << iter->id1 << " " << iter->id2 << " " << iter->tag << endl;
  }
  geofile << "EndConnect" << endl;
  geofile.close();
}
