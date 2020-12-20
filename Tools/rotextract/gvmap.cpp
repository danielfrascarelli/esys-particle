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

#include "gvmap.h"

#include <iterator>
#include <fstream>

using std::istream_iterator;
using std::back_inserter;
using std::ifstream;
using std::make_pair;

GVMap::GVMap(int tag)
{
  m_min_tag=tag;
}

void GVMap::read(const string& infilename, bool is_rot)
{
  ifstream headerfile(infilename.c_str()); 
  float dummy,xmax,ymax,zmax,xmin,ymin,zmin;
  string dummystring;
  vector<string> filenames; 
  int version;

  // clean multimap
  m_map.clear();

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
  // get main files
  for(vector<string>::iterator iter=filenames.begin();
      iter!=filenames.end();
      iter++){
    ifstream datafile(iter->c_str());
    vector<float> pdata;

    // get particles
    gdata gd;
    int npart;
    Vec3 oldpos;
    Vec3 initpos;
    Vec3 force;
    Vec3 angvel, circ_shift;
    double rad;
    int id;
    int tag;
    double q1,q2,q3,q4;
    datafile >> npart;
    for(int i=0;i<npart;i++){
      // read data
      if(is_rot){
	if(version<2){
	datafile >> gd.pos  >> rad >> id >> tag >> gd.mass >> initpos >> oldpos >> gd.vel >> force >> q1 >> q2 >> q3 >> q4 >> angvel;
      } else {
	datafile >> gd.pos  >> rad >> id >> tag >> gd.mass >> initpos >> oldpos >> gd.vel >> force >> circ_shift >> q1 >> q2 >> q3 >> q4 >> angvel;
      }
      } else {
	datafile >> gd.pos  >> rad >> id >> tag >> gd.mass >> initpos >> oldpos >> gd.vel >> force;
      }
      // put pos,vel,mass into map
      if(tag>=m_min_tag){
	m_map.insert(make_pair(tag,gd));
	// put tag into set
	m_tag_set.insert(tag);
      }
    }
  }
}

void GVMap::calc()
{
  for(set<int>::iterator iter=m_tag_set.begin();
      iter!=m_tag_set.end();
      iter++){
    Vec3 sum_pos=Vec3(0.0,0.0,0.0);
    Vec3 sum_vel=Vec3(0.0,0.0,0.0);
    double sum_mass=0;
    pair<mm_iter,mm_iter> iter_pair=m_map.equal_range(*iter);
    // 1st pass -> get grain center/velocity
    for(mm_iter it2=iter_pair.first;
	it2!=iter_pair.second;
	it2++){
      gdata g=it2->second;
      //      cout << g.pos << " - " << g.vel << " - " << g.mass << endl;
      // sum weigthed position
      sum_pos+=g.pos*g.mass;
      // sum weigthed velocities
      sum_vel+=g.vel*g.mass;
      // sum weight
      sum_mass+=g.mass;
    }
    // calc mass center of the grain
    Vec3 center=sum_pos/sum_mass;
    // calc vel. of mass center
    Vec3 grain_vel=sum_vel/sum_mass;
    cout << "mass, center, velocity of grain [ " << *iter << " ] : " << sum_mass << " ( " << center << " ) , ( " << grain_vel << " )" << endl;
    // 2nd pass -> get rotation
    Vec3 sum_angvel=Vec3(0.0,0.0,0.0);
    for(mm_iter it2=iter_pair.first;
	it2!=iter_pair.second;
	it2++){
      gdata g=it2->second;
      Vec3 r=g.pos-center;
      Vec3 angvel=cross(r,g.vel-grain_vel)/(r*r);
      cout << "[" << it2->first << "] " << angvel << " - " << g.mass << endl;
      sum_angvel+=angvel*g.mass;
    }
    Vec3 grain_angvel=sum_angvel/sum_mass;
    cout << "angular velocity of grain [ "  << *iter << " ] : " << grain_angvel << endl;
    m_vel[*iter]=grain_vel;
    m_angvel[*iter]=grain_angvel;
  }
}
