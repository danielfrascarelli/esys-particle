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

#include "fracframe.h"

//--- system includes ---
#include <map>
#include <fstream>
#include <iostream> 
#include <iterator>

using std::multimap;
using std::map;
using std::ifstream;
using std::cout;
using std::cerr;
using std::endl;
using std::istream_iterator;
using std::back_inserter;

//--- project includes ---
#include "Foundation/vec3.h"

FracFrame::FracFrame()
{}

bool FracFrame::cmp::operator()(const bdata& b1,const bdata& b2)
{
  bool res=false;

  if(b1.id1!=b2.id1){ 
    res=(b1.id1<b2.id1);
  } else {
    if (b1.id2!=b2.id2){ 
      res=(b1.id2<b2.id2);
    } else {
      res=(b1.tag<b2.tag);
    }
  }

  return res;
}

int FracFrame::get_version(const string& infilename)
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

vector<string> FracFrame::get_filenames(const string& infilename, int version)
{
  cout << "infilename : " << infilename << endl;
  ifstream headerfile(infilename.c_str()); 
  float dummy,xmax,ymax,zmax,xmin,ymin,zmin;
  vector<string> filenames;
  string dummystring;

  if(version==0){
    headerfile >> dummy >> dummy >> dummy;
    headerfile >> dummystring >> dummy;
  } else if ((version==1)||(version==2) || (version==3)) {
    headerfile >> dummystring >> dummy;
    headerfile >> dummy >> dummy >> dummy;
    std::getline(headerfile,dummystring);
    std::getline(headerfile,dummystring);
    std::cout << "dummystring:" << dummystring << std::endl;
  } else {
    cerr << "unknown checkpoint version " << version << endl;
  }
  // get bounding box
  headerfile >> dummystring;
  std::cout << "dummystring:" << dummystring << std::endl;
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
  read connection data from file

  \param infilename the filename
*/ 
void FracFrame::readFile(const string& infilename)
{
  int version=get_version(infilename.c_str());
  vector<string> filenames=get_filenames(infilename,version);


  // test 
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
    double rad;
    double mass;
    int id;
    int tag;
    datafile >> npart;
    if(version<2){
      for(int i=0;i<npart;i++){
	datafile >> pos  >> rad >> id >> tag >> mass >> initpos >> oldpos >> vel >> force;
	posmap[id]=pos;
	radmap[id]=rad;
      }
    } else {
      Vec3 circ_shift;
      for(int i=0;i<npart;i++){
	datafile >> pos  >> rad >> id >> tag >> mass >> initpos >> oldpos >> vel >> force >> circ_shift;
	posmap[id]=pos;
	radmap[id]=rad;
      }
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
	  bdata bond;
	  
	  datafile >> bond.id1 >> bond.id2 >> bond.tag;
	  m_bonds.insert(bond);
	}
      }
    } else { // pre - V1 snapshot -> assume bondend pair IG
      datafile >> nbond;
      for(int i=0;i<nbond;i++){
	bdata bond;
	
	datafile >> bond.id1 >> bond.id2 >> bond.tag;
	m_bonds.insert(bond);
      }
    }



    datafile.close();
  }
}

/*!
  read connection data involving particles with a particular tag from file

  \param infilename the filename
  \param rtag the particle tag
*/ 
void FracFrame::readFileTagged(const string& infilename, int rtag)
{
  int version=get_version(infilename.c_str());
  vector<string> filenames=get_filenames(infilename,version);


  // test 
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
    double rad;
    double mass;
    int id;
    int tag;
    datafile >> npart;
    if(version<2){
      for(int i=0;i<npart;i++){
		datafile >> pos  >> rad >> id >> tag >> mass >> initpos >> oldpos >> vel >> force;
		posmap[id]=pos;
		radmap[id]=rad;
		tagmap[id]=tag;
	  }
    } else {
      Vec3 circ_shift;
      for(int i=0;i<npart;i++){
		datafile >> pos  >> rad >> id >> tag >> mass >> initpos >> oldpos >> vel >> force >> circ_shift;
		posmap[id]=pos;
		radmap[id]=rad;
		tagmap[id]=tag;
	  }
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
	  bdata bond;
	  
	  datafile >> bond.id1 >> bond.id2 >> bond.tag;
	  m_bonds.insert(bond);
	}
      }
    } else { // pre - V1 snapshot -> assume bondend pair IG
      datafile >> nbond;
      for(int i=0;i<nbond;i++){
	bdata bond;
	
	datafile >> bond.id1 >> bond.id2 >> bond.tag;
	m_bonds.insert(bond);
      }
    }



    datafile.close();
  }
}

/*!
  read connection data from a file (rotational particles)

  \param infilename the filename
  \param init_pos if true, write fractures at _initial_ particle positions
*/
void FracFrame::readFileRot(const string& infilename,bool init_pos)
{
  int version=get_version(infilename.c_str());
  vector<string> filenames=get_filenames(infilename,version);


  // test 
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
		if(init_pos) {
			posmap[id]=initpos;
		} else {
			posmap[id]=pos;
		}
		radmap[id]=rad;
		tagmap[id]=tag;
	  }
    } else {
      Vec3 circ_shift;
      for(int i=0;i<npart;i++){
		datafile >> pos  >> rad >> id >> tag >> mass >> initpos >> oldpos >> vel >> force >> circ_shift >> q1 >> q2 >> q3 >> q4 >> angvel;
		if(init_pos) {
		  posmap[id]=initpos;
		} else {
		  posmap[id]=pos;
		}
		radmap[id]=rad;
		tagmap[id]=tag;
     }
    }
      
    cout << "nr. of particles: " << posmap.size() << endl;

    // get bonds
    int ngrp;
    int nbond;
    string type;

    datafile >> ngrp;
    if(version>0){
      for (int j=0;j<ngrp; j++){
	datafile >> type;
	if((type=="Bonded") || (type=="RotBonded")){
	  datafile >> nbond;
	  cout << "nbond: " << nbond << endl;
	  for(int i=0;i<nbond;i++){
	    bdata bond;
	    
	    datafile >> bond.id1 >> bond.id2 >> bond.tag;
	    m_bonds.insert(bond);
	  }
	}
      }
    } else { // pre - V1 snapshot -> assume bondend pair IG
      datafile >> nbond;
      for(int i=0;i<nbond;i++){
		bdata bond;
			
		datafile >> bond.id1 >> bond.id2 >> bond.tag;
		m_bonds.insert(bond);
      }
    }

    cout << "nr. of bonds is:" << m_bonds.size() << endl;

    datafile.close();
  }
}

/*!
  read connection data involving particles with a particular tag from a file (rotational particles)

  \param infilename the filename
  \param rtag the particle tag
*/
void FracFrame::readFileRotTagged(const string& infilename,int rtag)
{
  int version=get_version(infilename.c_str());
  vector<string> filenames=get_filenames(infilename,version);
  set<int> idset; // set of ids of particles with the correct tag 

  // test 
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
	if(tag==rtag){
	  posmap[id]=pos;
	  radmap[id]=rad;
	  tagmap[id]=tag;
	  idset.insert(id);
	}
      }
    } else {
      Vec3 circ_shift;
      for(int i=0;i<npart;i++){
	datafile >> pos  >> rad >> id >> tag >> mass >> initpos >> oldpos >> vel >> force >> circ_shift >> q1 >> q2 >> q3 >> q4 >> angvel;
	if(tag==rtag){
	  posmap[id]=pos;
	  radmap[id]=rad; 
	  tagmap[id]=tag;
	  idset.insert(id);
	}
      }
    }
      
    cout << "nr. of particles: " << posmap.size() << endl;

    // get bonds
    int ngrp;
    int nbond;
    string type;

    datafile >> ngrp;
    if(version>0){
      datafile >> type;
      if((type=="Bonded") || (type=="RotBonded")){
	datafile >> nbond;
	cout << "nbond: " << nbond << endl;
	for(int i=0;i<nbond;i++){
	  bdata bond;
	    
	  datafile >> bond.id1 >> bond.id2 >> bond.tag;
	  if((idset.find(bond.id1)!=idset.end()) || (idset.find(bond.id2)!=idset.end())){
	    m_bonds.insert(bond);	  
	  }
	}
      }
    } else { // pre - V1 snapshot -> assume bondend pair IG
      datafile >> nbond;
      for(int i=0;i<nbond;i++){
	bdata bond;
	    
	datafile >> bond.id1 >> bond.id2 >> bond.tag;
	if((idset.find(bond.id1)!=idset.end()) || (idset.find(bond.id2)!=idset.end())){
	  m_bonds.insert(bond);	  
	}
      }
    }

    cout << "nr. of bonds is:" << m_bonds.size() << endl;

    datafile.close();
  }
}

/*!
  get the fracture as difference in the bonds between two frames

  \param F2 the newer frame
*/
vector<FracFrame::fdata> FracFrame::getFrac(FracFrame& F2)
{
  vector<fdata> res;
  
  // debug output
  std::cout << "nr. of particles old: " << radmap.size() << " bonds: " << m_bonds.size() << std::endl;
  std::cout << "nr. of particles new: " << F2.radmap.size() << " bonds: " << F2.m_bonds.size() << std::endl;


  // find bonds in this but not in F2
  for(set<bdata>::const_iterator iter=m_bonds.begin();
      iter!=m_bonds.end();
      iter++){
    if((F2.m_bonds).find(*iter)==(F2.m_bonds).end()){
      fdata nfd;

      cout << "bond [" << iter->id1 << "-" << iter->id2 << "] broken" << endl; 
      // check if both particles are in new frame
      map<int,float>::iterator op1=F2.radmap.find(iter->id1);
      if(op1==F2.radmap.end()) std::cout << "particle " << iter->id1 << " gone" << std::endl;
      map<int,float>::iterator op2=F2.radmap.find(iter->id2);
      if(op2==F2.radmap.end()) std::cout << "particle " << iter->id2 << " gone" << std::endl;
      
      if((op1!=F2.radmap.end())&&(op2!=F2.radmap.end())){
		// particle posn
		Vec3 pos1=F2.posmap[iter->id1];
		Vec3 pos2=F2.posmap[iter->id2];
		// particle rad
		double r1=F2.radmap[iter->id1];
		double r2=F2.radmap[iter->id2];
		// contact normal
		nfd.normal=(pos2-pos1).unit();
		// contact pos.
		double compr=(pos2-pos1).norm()/(r1+r2);
		nfd.pos=pos1+r1*compr*nfd.normal;
		// size
		nfd.size=0.5*(r1+r2);
		// distance - can be used to filter out "circular" breaks
		nfd.dist=(pos2-pos1).norm();
		nfd.id1=iter->id1;
		nfd.id2=iter->id2;
		nfd.tag=iter->tag;
		if (tagmap[iter->id1] <= tagmap[iter->id2]){
			nfd.ptag1=tagmap[iter->id1];
			nfd.ptag2=tagmap[iter->id2];
		} else {	
			nfd.ptag1=tagmap[iter->id2];
			nfd.ptag2=tagmap[iter->id1];
		}
		res.push_back(nfd);
      }
    }
  }

  return res;
}
