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

//--- project include ---
#include "../../Foundation/vec3.h"
#include "rextract.h"

//--- STL includes ---
#include <sstream>
#include <iterator>
#include <fstream>
#include <vector>
#include <cmath>

using std::ostringstream;
using std::istream_iterator;
using std::back_inserter;
using std::ifstream;
using std::make_pair;
using std::vector;
using std::sqrt;

Rextract::Rextract(const string& ifn,const string& ofn,int tag)
{
  m_infilebase=ifn;
  m_outfilename=ofn;
  m_tag=tag;
  m_count=0;
  m_sumpart=0;
  m_initialized=false;
}


void Rextract::read_frame(int t)
{
  // make filename
  ostringstream infilename;
  infilename << m_infilebase << "_t=" << t << "_0.txt";
  std::cout << "header file: " << infilename.str() << std::endl; 
  ifstream headerfile((infilename.str()).c_str()); 

  float dummy,xmax,ymax,zmax,xmin,ymin,zmin;
  string dummystring;
  vector<string> filenames; 
  int version;

  // read headerfile
  headerfile >> dummystring;
  if(dummystring=="V"){ // if V -> new version 
    headerfile >> version ;
    cout << "version : " << version << endl;
    headerfile >> dummy >> dummy >> dummy; 
  } else {
    cout << "pre- V.1 version" << endl;
    version=0;
    headerfile >> dummy >> dummy;
  }
  headerfile >> dummystring >> dummy; // geometry version

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

  // if first frame, get particle numbers
  if(!m_initialized){
    for(vector<string>::iterator iter=filenames.begin();
	iter!=filenames.end();
	iter++){
      int npart;
      ifstream datafile(iter->c_str());
      datafile >> npart;
      m_sumpart+=npart;
      datafile.close();
    }
    m_angvel_rms=vector<double>(m_sumpart);
    m_rad=vector<double>(m_sumpart);
    m_initialized=true;
  }

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
  // get main files
  for(vector<string>::iterator iter=filenames.begin();
      iter!=filenames.end();
      iter++){
    ifstream datafile(iter->c_str());
   
    datafile >> npart;
    std::cout << "file: " << *iter << " particles: " << npart << std::endl;

    for(int i=0;i<npart;i++){
      // read data
      if(version<2){
	datafile >> pos  >> rad >> id >> tag >> mass >> initpos >> oldpos >> vel >> force >> q1 >> q2 >> q3 >> q4 >> angvel;
      } else {
	datafile >> pos  >> rad >> id >> tag >> mass >> initpos >> oldpos >> vel >> force >> circ_shift >> q1 >> q2 >> q3 >> q4 >> angvel;
      }
      // if correct tag -> add angvel^2 
      if(tag==m_tag){
	m_angvel_rms[id]+=angvel*angvel;
	if(m_count==0){
	  m_rad[id]=rad;
	  m_tag_set.insert(id);
	}
      }
    }
    datafile.close();
  }
  m_count++;
}

void Rextract::write_data()
{
  ofstream outfile(m_outfilename.c_str());
  for(set<int>::iterator iter=m_tag_set.begin();
      iter!=m_tag_set.end();
      iter++){
    outfile << m_rad[*iter] << "  " << sqrt(m_angvel_rms[*iter]/double(m_count)) << endl;
  }
  outfile.close();
}

void Rextract::write_data_bin(double x0,double dx,int nb)
{
  vector<double> r_2(nb,0.0);
  vector<double> r_m(nb,0.0);
  vector<double> omega_2(nb,0.0);
  vector<double> omega_m(nb,0.0);
  vector<int> count(nb,0);

  ofstream outfile(m_outfilename.c_str());
  for(set<int>::iterator iter=m_tag_set.begin();
      iter!=m_tag_set.end();
      iter++){
    double r=m_rad[*iter];
    double omega=sqrt(m_angvel_rms[*iter]/double(m_count));
    int idx=int(floor((r-x0)/dx));
    if((idx>=0) && (idx<nb)){
      r_m[idx]+=r;
      omega_m[idx]+=omega;
      count[idx]++;
    }
  }
  for(int i=0;i<nb;i++){
    outfile << r_m[i]/double(count[i]) << " " << omega_m[i]/double(count[i]) << endl;
  }
  outfile.close();
}
