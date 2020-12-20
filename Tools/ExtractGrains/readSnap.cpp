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

#include "readSnap.h"

#include <iterator>
#include <fstream>
#include <set>

using std::istream_iterator;
using std::back_inserter;
using std::ifstream;
using std::make_pair;
using std::set;

//--- Project includes ---
#include "Foundation/vec3.h"

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

vector<string> get_filenames(const string& infilename, int version)
{
  cout << "infilename : " << infilename << endl;
  ifstream headerfile(infilename.c_str()); 
  float dummy,xmax,ymax,zmax,xmin,ymin,zmin;
  vector<string> filenames;
  string dummystring;

  if(version==0){
    headerfile >> dummy >> dummy >> dummy;
    headerfile >> dummystring >> dummy;
  } else if ((version==1) || (version==2) || (version==3)){
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

  cout << "nr. of filenames: " << filenames.size() << endl;
  return filenames;
}

/*!
  read snapshot and create a graph from the data

  \param infilename the name of the header file of the snapshot, i.e. *_0.txt
  \param mintag the smallest particle tag considered
  \param btag if btag!= -1, consider only bonds tagged with btag
*/
Graph readSnap(const string& infilename,int mintag,int maxtag,int btag)
{
  ifstream headerfile(infilename.c_str());
  Graph g;

  //float dummy=0.0,xmax=0.0,ymax=0.0,zmax=0.0,xmin=0.0,ymin=0.0,zmin=0.0;
  string dummystring;
  set<int> idset;
  set<pair<int,int> > edgeset;

  int version=get_version(infilename);
  vector<string> filenames=get_filenames(infilename,version);

  // get main files
  for(vector<string>::iterator iter=filenames.begin();
      iter!=filenames.end();
      iter++){
    cout << *iter << endl;
    ifstream datafile(iter->c_str());

    // get particles
    int npart;
    Vec3 pos;
    Vec3 oldpos;
    Vec3 initpos;
    Vec3 force;
    Vec3 vel;
    Vec3 angvel;
    Vec3 circ_shift;
    double rad;
    double mass;
    double q1,q2,q3,q4;
    int id;
    int tag;
    datafile >> npart;
    std::cerr << "parts :" << npart << std::endl;
    if(version<2){
      for(int i=0;i<npart;i++){
	// read data
	datafile >> pos  >> rad >> id >> tag >> mass >> initpos >> oldpos >> vel >> force >> q1 >> q2 >> q3 >> q4 >> angvel;
	// if correct tag -> add id to set
	if((tag>=mintag) && (tag<=maxtag)) {
	  idset.insert(id);
	  g.setVertexData(id,pdata(tag,mass,pos,rad,vel,angvel));
	}
      }
    } else {
      for(int i=0;i<npart;i++){
	// read data
	datafile >> pos  >> rad >> id >> tag >> mass >> initpos >> oldpos >> vel >> force >> circ_shift >> q1 >> q2 >> q3 >> q4 >> angvel;
	// if correct tag -> add id to set
	if((tag>=mintag) && (tag <= maxtag)){
	  idset.insert(id);
	  g.setVertexData(id,pdata(tag,mass,pos,rad,vel,angvel));
	}
      }
    }
    
    // get bonds

    int ngroup;
    int nbond;
    int p1,p2,bondtag;
    string type;

    datafile >> ngroup;
    std::cerr << "groups :" << ngroup << std::endl;
    if(version>0){
      for(int i=0;i<ngroup;i++){
	datafile >> type;
	if((type=="Bonded") || (type=="RotBonded")){
	  datafile >> nbond;
	  std::cerr << "bonds :" << nbond << std::endl;
	  for(int j=0;j<nbond;j++){
	    datafile >> p1 >> p2 >> bondtag;
	    if((btag==-1) || (bondtag==btag)){ // if bond tag fits or is ignored
	      if((idset.find(p1)!=idset.end()) && (idset.find(p2)!=idset.end())){
		edgeset.insert(make_pair(p1,p2));
	      } else {
		// std::cerr << "bond with missing particle: " << p1 << " - " << p2 << std::endl; 
	      }
	    }
	  }
	}
      }
    } else  { // pre - V1 snapshot -> assume bondend pair IG
      for(int i=0;i<ngroup;i++){
	datafile >> nbond;
	for(int j=0;j<nbond;j++){
	  datafile >> p1 >> p2 >> btag;
	  if((btag==-1) || (bondtag==btag)){ // if bond tag fits or is ignored
	    if((idset.find(p1)!=idset.end()) && (idset.find(p2)!=idset.end())){
	      edgeset.insert(make_pair(p1,p2));
	    }
	  }
	}
      }
    } 
  }
   
  for(set<pair<int,int> >::iterator iter=edgeset.begin();
      iter!=edgeset.end();
      iter++){
    g.insert(*iter);
  }
  cout << "size of idset: " << idset.size() << "  size of edgeset : " << edgeset.size() << endl;
  cout << "nr. of vertices and edges: " << g.numV() << " , " << g.numE() << endl;
  g.removeDoubles();
  cout << "nr. of vertices and edges after removal of doubles: " << g.numV() << " , " << g.numE() << endl;
  cout << endl;
  return g;
}

/*!
  read snapshot into existing graph

  \param infilename the name of the header file of the snapshot, i.e. *_0.txt
  \param mintag the smallest particle tag considered
  \param g the existing graph
  \param btag if btag!= -1, consider only bonds tagged with btag
*/
void readSnap(const string& infilename,int mintag, int maxtag, Graph& g, int btag)
{
  ifstream headerfile(infilename.c_str());
  std::cerr << "bond tag: " << btag << std::endl;
  //float dummy=0.0,xmax=0.0,ymax=0.0,zmax=0.0,xmin=0.0,ymin=0.0,zmin=0.0;
  string dummystring;
  set<int> idset;
  set<pair<int,int> > edgeset;

  int version=get_version(infilename);
  vector<string> filenames=get_filenames(infilename,version);

  // get main files
  for(vector<string>::iterator iter=filenames.begin();
      iter!=filenames.end();
      iter++){
    cout << *iter << endl;
    ifstream datafile(iter->c_str());

    // get particles
    int npart;
    Vec3 pos;
    Vec3 oldpos;
    Vec3 initpos;
    Vec3 force;
    Vec3 vel;
    Vec3 angvel;
    Vec3 circ_shift;
    double rad;
    double mass;
    double q1,q2,q3,q4;
    int id;
    int tag;
    datafile >> npart;
    std::cerr << "parts :" << npart << std::endl;
    if(version<2){
      for(int i=0;i<npart;i++){
	// read data
	datafile >> pos  >> rad >> id >> tag >> mass >> initpos >> oldpos >> vel >> force >> q1 >> q2 >> q3 >> q4 >> angvel;
	// if correct tag -> add id to set
	if((tag>=mintag) && (tag <= maxtag)) {
	  idset.insert(id);
	  g.setVertexData(id,pdata(tag,mass,pos,rad,vel,angvel));
	}
      }
    } else {
      for(int i=0;i<npart;i++){
	// read data
	datafile >> pos  >> rad >> id >> tag >> mass >> initpos >> oldpos >> vel >> force >> circ_shift >> q1 >> q2 >> q3 >> q4 >> angvel;
	// if correct tag -> add id to set
	if((tag>=mintag) && (tag <= maxtag)){
	  idset.insert(id);
	  g.setVertexData(id,pdata(tag,mass,pos,rad,vel,angvel));
	}
      }
    }
    
    // get bonds

    int ngroup;
    int nbond;
    int p1,p2,bondtag;
    string type;

    datafile >> ngroup;
    std::cerr << "groups :" << ngroup << std::endl;
    if(version>0){
      for(int i=0;i<ngroup;i++){
	datafile >> type;
	if((type=="Bonded") || (type=="RotBonded")){
	  datafile >> nbond;
	  std::cerr << "bonds :" << nbond << std::endl;
	  for(int j=0;j<nbond;j++){
	    datafile >> p1 >> p2 >> bondtag;
	    if((btag==-1) || (bondtag==btag)){ // if bond tag fits or is ignored
	      edgeset.insert(make_pair(p1,p2));
	    }
	  }
	}
      }
    } else  { // pre - V1 snapshot -> assume bondend pair IG
      for(int i=0;i<ngroup;i++){
	datafile >> nbond;
	for(int j=0;j<nbond;j++){
	  datafile >> p1 >> p2 >> bondtag;
	  if((btag==-1) || (bondtag==btag)){ // if bond tag fits or is ignored
	    if((idset.find(p1)!=idset.end()) && (idset.find(p2)!=idset.end())){
	      edgeset.insert(make_pair(p1,p2));
	    }
	  }
	}
      }
    } 
  }
   
  for(set<pair<int,int> >::iterator iter=edgeset.begin();
      iter!=edgeset.end();
      iter++){
    if((idset.find(iter->first)!=idset.end()) && (idset.find(iter->second)!=idset.end())){
      g.insert(*iter);
    }
  }
  cout << "size of idset: " << idset.size() << "  size of edgeset : " << edgeset.size() << endl;
  cout << "nr. of vertices and edges: " << g.numV() << " , " << g.numE() << endl;
  g.removeDoubles();
  cout << "nr. of vertices and edges after removal of doubles: " << g.numV() << " , " << g.numE() << endl;
  cout << endl;
}

/*!
  Read geometry file and create a graph from the data

  \param infilename the name of the geometry file 
  \param mintag the smallest particle tag considered
  \param btag if btag!= -1, consider only bonds tagged with btag
*/
Graph readGeo(const string& infilename,int mintag,int maxtag,int btag)
{
  Graph g;

  return g;
}


/*!
  Read geometry file into existing graph

  \param infilename the name of the geometry file 
  \param mintag the smallest particle tag considered
  \param g the existing graph
  \param btag if btag!= -1, consider only bonds tagged with btag
*/
void readGeo(const string& infilename,int mintag,int maxtag,Graph& g,int btag)
{
  string dummy;
  Vec3 pos;
  double rad;
  int nparts,nbonds;
  int id;
  int tag;

  set<int> idset;
  set<pair<int,int> > edgeset;

  // open file
  ifstream ifile(infilename.c_str());

  // read and forget 6 lines
  ifile >> dummy >> dummy;
  ifile >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy ;
  ifile >> dummy >> dummy >> dummy >> dummy ;
  ifile >> dummy >> dummy;
  ifile >> dummy;
  ifile >> dummy;
  // number of particles
  ifile >> nparts;

  // read particles
  std::cout << "nr. of particles " << nparts << std::endl;
  for(int i=0;i<nparts;i++){
    ifile >> pos >> rad >> id >> tag;
    // if correct tag -> add id to set
    if((tag>=mintag) && (tag <=maxtag)){
      idset.insert(id);
      double mass=(4.0/3.0)*3.14159*rad*rad*rad;
      g.setVertexData(id,pdata(tag,mass,pos,rad,Vec3(0.0,0.0,0.0),Vec3(0.0,0.0,0.0)));
    }
    
  }
  // read bonds
  ifile >> dummy >> dummy;
  ifile >> nbonds;
  std::cout << "nr. of bonds " << nbonds << std::endl;
  int p1,p2,bondtag;
  for(int j=0;j<nbonds;j++){
    ifile >> p1 >> p2 >> bondtag;
    if((btag==-1) || (bondtag==btag)){ // if bond tag fits or is ignored
      if((idset.find(p1)!=idset.end()) && (idset.find(p2)!=idset.end())){
	edgeset.insert(make_pair(p1,p2));
      } else { 
	std::cout << "no insert: " << p1 << " [" << (idset.find(p1)!=idset.end()) << " " << p2 <<   " [" << (idset.find(p2)!=idset.end()) << "]" << std::endl;
      }
    }
  }
  for(set<pair<int,int> >::iterator iter=edgeset.begin();
      iter!=edgeset.end();
      iter++){
    if((idset.find(iter->first)!=idset.end()) && (idset.find(iter->second)!=idset.end())){
      g.insert(*iter);
    }
  }
  cout << "size of idset: " << idset.size() << "  size of edgeset : " << edgeset.size() << endl;
  cout << "nr. of vertices and edges: " << g.numV() << " , " << g.numE() << endl;
  g.removeDoubles();
  cout << "nr. of vertices and edges after removal of doubles: " << g.numV() << " , " << g.numE() << endl;
  cout << endl;
  // close file
  ifile.close();
  
}
