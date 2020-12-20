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

// --- TML includes ---
#include "tml/comm/comm.h"

// --- Project includes ---
#include "ParticleFieldMaster.h"
#include "field_const.h"

//--- IO includes ---
#include <iostream>
#include <fstream>

//--- STL inculdes ---
#include <string>
#include <map>

using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::multimap;

//==== SCALAR PFM ===

/*!
  Constructor. Set up the Master and broadcast parameters to the slaves.

  \param comm the communicator
  \param fieldname the name of the field to be saved
  \param filename the name of the file to be saved into or the base for the generation of the filenames if the saving format requires multiple files
  \param savetype the way to save data, currently supported are DX,SUM
  \param t0 the first timestep to be saved
  \param tend the last timestep to be saved
  \param dt the interval between timesteps to be saving
*/ 
ScalarParticleFieldMaster::ScalarParticleFieldMaster(TML_Comm* comm,const string& fieldname,const string& filename,const string& savetype,int t0,int tend,int dt)
  :AFieldMaster(comm,fieldname,filename,savetype,t0,tend,dt)
{
  m_comm->broadcast_cont(fieldname);
  m_comm->broadcast(m_id);
  int is_tagged=0;
  m_comm->broadcast(is_tagged);
}

/*!
  Constructor. Set up the Master and broadcast parameters to the slaves.

  \param comm the communicator
  \param fieldname the name of the field to be saved
  \param filename the name of the file to be saved into or the base for the generation of the filenames if the saving format requires multiple files
  \param savetype the way to save data, currently supported are DX,SUM
  \param t0 the first timestep to be saved
  \param tend the last timestep to be saved
  \param dt the interval between timesteps to be saving
  \param tag the tag of the particles to be saved
  \param mask the mask to be applied to the tag
*/ 
ScalarParticleFieldMaster::ScalarParticleFieldMaster(TML_Comm* comm,const string& fieldname,const string& filename,const string& savetype,int t0,int tend,int dt,int tag,int mask):AFieldMaster(comm,fieldname,filename,savetype,t0,tend,dt)
{
  m_comm->broadcast_cont(fieldname);
  m_comm->broadcast(m_id);
  int is_tagged=1;
  m_comm->broadcast(is_tagged);
  m_comm->broadcast(tag);
  m_comm->broadcast(mask);
}



void ScalarParticleFieldMaster::collect()
{
  // send field id to slave
  m_comm->broadcast(m_id);

  if((m_write_type==WRITE_TYPE_SUM)||(m_write_type==WRITE_TYPE_MAX)) {
    collectSum(); 
  } else {
    collectFull();
  }

}

/*! 
  collect the full set of data, i.e. value, and posn,size of particles
*/
void ScalarParticleFieldMaster::collectFull()
{ 
  multimap<int,pair<int,double> > temp_mm;
  multimap<int,pair<int,Vec3> > temp_pos_mm;
  multimap<int,pair<int,double> > temp_rad_mm;
  
  // send type of collect to slaves
  int coll_type=1;
  m_comm->broadcast(coll_type);

  // get data from slaves
  m_comm->gather(temp_mm);
  m_comm->gather(temp_pos_mm);
  m_comm->gather(temp_rad_mm);

  // collate receive data from multimap into single map  
  for(multimap<int,pair<int,double> >::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    m_save_map.insert(iter->second);
  }
  for(multimap<int,pair<int,Vec3> >::iterator iter=temp_pos_mm.begin();
      iter!=temp_pos_mm.end();
      iter++){
    m_pos_map.insert(iter->second);
  }
  for(multimap<int,pair<int,double> >::iterator iter=temp_rad_mm.begin();
      iter!=temp_rad_mm.end();
      iter++){
    m_rad_map.insert(iter->second);
  }
}


/*!
  collect sum of values only 
*/
void ScalarParticleFieldMaster::collectSum()
{
  multimap<int,double> temp_mm;

  // send type of collect to slaves
  m_comm->broadcast(m_write_type);

  // get data from slaves
  m_comm->gather(temp_mm);

  // collate receive data from multimap into single map  
  for(multimap<int,double>::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    m_save_map.insert(*iter);
  } 
}


/*!
  write data out as OpenDX compatible file

  \todo desciption
*/
void ScalarParticleFieldMaster::writeAsDX()
{
  //generate filename
  string fn=makeFilename();
  // write header 
  ofstream out_file(fn.c_str());
  out_file << "points = " << m_save_map.size() << endl;
  out_file << "format = ascii" << endl;
  out_file << "dependency = positions, positions" << endl;
  out_file << "interleaving = field" << endl;
  out_file << "field = locations, " << m_field_name << endl;
  out_file << "structure = 3-vector, scalar" << endl;
  out_file << "type = float, float  " << endl;
  out_file << "header =  marker \"Start\\n\"" << endl;
  out_file << endl << "end" << endl;
  out_file << "Start" << endl;
  
  // write data
  for(map<int,double>::iterator iter=m_save_map.begin();
      iter!=m_save_map.end();
      iter++){
    out_file << m_pos_map[iter->first]  << " " << iter->second << endl;
  }
  out_file.close();
  m_save_map.erase(m_save_map.begin(),m_save_map.end());
  m_pos_map.erase(m_pos_map.begin(),m_pos_map.end());
  m_rad_map.erase(m_rad_map.begin(),m_rad_map.end());
}

/*!
  write data out as Povray(>=3.0) compatible file
  
  \warning not fully impl.
*/
void ScalarParticleFieldMaster::writeAsPOV()
{
  //generate filename
  string fn=makeFilename();
  // open file
  ofstream out_file(fn.c_str());

  // work out object size
  Vec3 vmin,vmax;
  for(map<int,Vec3>::iterator iter=m_pos_map.begin();
      iter!=m_pos_map.end();
      iter++){
    vmin=cmin(vmin,iter->second);
    vmax=cmax(vmax,iter->second);
  }
  // fit camera position
  //Vec3 lookat=0.5*(vmin+vmax);
  
  // write camera etc.
  out_file << "#include \"colors.inc\"" << endl;
  out_file << "background { color Black }" << endl;
  out_file << " camera {" << endl;
  out_file << "location <0, 2, -3>" << endl;
  out_file << "look_at  <0, 1,  2>" << endl;
  out_file << "}" << endl<< endl;

  out_file << "light_source { <2, 4, -3> color White}" << endl;

  // write particles
  for(map<int,double>::iterator iter=m_save_map.begin();
      iter!=m_save_map.end();
      iter++){
    out_file << "sphere { <" << m_pos_map[iter->first].X() << " , "
	     << m_pos_map[iter->first].Y() << " , " 
	     << m_pos_map[iter->first].Z() << " > , "
	     << m_rad_map[iter->first] << endl;
    out_file << " texture { rgb < 1.0,0.5,0.5 > } }" << endl; //<< iter->second << endl;
  }
  // close file 
  out_file.close();
}

/*!
  write data out as SILO file (if supported)
*/
void ScalarParticleFieldMaster::writeAsSILO()
{
#if HAVE_SILO
  bool append;
  DBfile* dbfile = openSiloFile(append);

  if (!dbfile) {
    console.Error() << "writeAsSILO not supported";
    return;
  }

  // prepare data
  std::vector<float> x, y, z, radius, data;
  x.clear(); y.clear(); z.clear(); radius.clear(); data.clear();

  for (map<int,double>::iterator iter=m_save_map.begin();
       iter != m_save_map.end();
       iter++) {
    // the mesh is only needed if we are not appending
    if (!append) {
      const Vec3& v = m_pos_map[iter->first];
      x.push_back(v.X());
      y.push_back(v.Y());
      z.push_back(v.Z());
    }
    radius.push_back(m_rad_map[iter->first]);
    data.push_back(iter->second);
  }

  // write mesh first if needed
  if (!append) {
    float* coords[] = { &x[0], &y[0], &z[0] };
    DBPutPointmesh(dbfile, "mesh", 3, coords, m_save_map.size(),
                   DB_FLOAT, NULL);
  }

  // and now the fields
  DBPutPointvar1(dbfile, "radius", "mesh", &radius[0], radius.size(),
                 DB_FLOAT, NULL);
  DBPutPointvar1(dbfile, m_field_name.c_str(), "mesh", &data[0],
                 data.size(), DB_FLOAT, NULL);
  // close file 
  DBClose(dbfile);
  m_save_map.erase(m_save_map.begin(),m_save_map.end());
  m_pos_map.erase(m_pos_map.begin(),m_pos_map.end());
  m_rad_map.erase(m_rad_map.begin(),m_rad_map.end());
#else
  console.Error() << "writeAsSILO not supported\n";
#endif
}

/*!
  sum data and write them out into a single continuous file

  \warning n
*/
void ScalarParticleFieldMaster::writeAsSUM()
{
  // sum data
  double sum_data=0.0;
  for(map<int,double>::iterator iter=m_save_map.begin();
      iter!=m_save_map.end();
      iter++){
    sum_data+=iter->second;
  }
  // open file for appending
  ofstream out_file(m_file_name.c_str(),ios::app);
  // write data
  out_file << sum_data << endl;
  // close file
  out_file.close();
  m_save_map.erase(m_save_map.begin(),m_save_map.end());
  m_pos_map.erase(m_pos_map.begin(),m_pos_map.end());
  m_rad_map.erase(m_rad_map.begin(),m_rad_map.end());

}

/*!
  find max datum and write it out into a single continuous file
*/
void ScalarParticleFieldMaster::writeAsMAX()
{
  // get max data
  double max_data=(m_save_map.begin())->second;
  for(map<int,double>::iterator iter=m_save_map.begin();
      iter!=m_save_map.end();
      iter++){
   	max_data=(iter->second > max_data) ? iter->second : max_data;
  }
  // open file for appending
  ofstream out_file(m_file_name.c_str(),ios::app);
  // write data
  out_file << max_data << endl;
  // close file
  out_file.close();
  m_save_map.erase(m_save_map.begin(),m_save_map.end());
  m_pos_map.erase(m_pos_map.begin(),m_pos_map.end());
  m_rad_map.erase(m_rad_map.begin(),m_rad_map.end());

}

/*!
  write data as a raw series of values, one row of values per time step,
  all timesteps into the same file
*/
void ScalarParticleFieldMaster::writeAsRAW_SERIES()
{  
  // open file
  ofstream out_file(m_file_name.c_str(),ios::app);
  
  // write data
  for(map<int,double>::iterator iter=m_save_map.begin();
      iter!=m_save_map.end();
      iter++){
    out_file << iter->second << "  ";
  }
  out_file << endl; 
  out_file.close();
  m_save_map.erase(m_save_map.begin(),m_save_map.end());
  m_pos_map.erase(m_pos_map.begin(),m_pos_map.end());
  m_rad_map.erase(m_rad_map.begin(),m_rad_map.end());
}

/*!
  write data as raw id, position, radius, value, one time step per file
*/
void ScalarParticleFieldMaster::writeAsRawWithPosID()
{
  //generate filename
  string fn=makeFilename();

  // open file
  ofstream out_file(fn.c_str());

  // write data
  for(map<int,double>::iterator iter=m_save_map.begin();
      iter!=m_save_map.end();
      iter++){
    out_file << iter->first << " " << m_pos_map[iter->first] << " " << m_rad_map[iter->first] << " "<< iter->second << endl;
  }

  // close file
  out_file.close();

  // clean up
  m_save_map.erase(m_save_map.begin(),m_save_map.end());
  m_pos_map.erase(m_pos_map.begin(),m_pos_map.end());
  m_rad_map.erase(m_rad_map.begin(),m_rad_map.end());
}

// === VECTOR PFM ===

/*!
  Constructor. Set up the Master and broadcast parameters to the slaves.

  \param comm the communicator
  \param fieldname the name of the field to be saved
  \param filename the name of the file to be saved into or the base for the generation of the filenames if the saving format requires multiple files
  \param savetype the way to save data, currently supported are DX,SUM
  \param t0 the first timestep to be saved
  \param tend the last timestep to be saved
  \param dt the interval between timesteps to be saving
*/ 
VectorParticleFieldMaster::VectorParticleFieldMaster(TML_Comm* comm,const string& fieldname,const string& filename,const string& savetype,int t0,int tend,int dt)
  :AFieldMaster(comm,fieldname,filename,savetype,t0,tend,dt)
{
  m_comm->broadcast_cont(fieldname);
  m_comm->broadcast(m_id);
  int is_tagged=0;
  m_comm->broadcast(is_tagged);
}

/*!
  Constructor. Set up the Master and broadcast parameters to the slaves.

  \param comm the communicator
  \param fieldname the name of the field to be saved
  \param filename the name of the file to be saved into or the base for the generation of the filenames if the saving format requires multiple files
  \param savetype the way to save data, currently supported are DX,SUM
  \param t0 the first timestep to be saved
  \param tend the last timestep to be saved
  \param dt the interval between timesteps to be saving
  \param tag the tag of the particles to be saved
  \param mask the mask to be applied to the tag
*/ 
VectorParticleFieldMaster::VectorParticleFieldMaster(TML_Comm* comm,const string& fieldname,const string& filename,const string& savetype,int t0,int tend,int dt,int tag,int mask)
  :AFieldMaster(comm,fieldname,filename,savetype,t0,tend,dt)
{
  m_comm->broadcast_cont(fieldname);
  m_comm->broadcast(m_id);
  int is_tagged=1;
  m_comm->broadcast(is_tagged);
  m_comm->broadcast(tag);
  m_comm->broadcast(mask);
}

void VectorParticleFieldMaster::collect()
{
  multimap<int,pair<int,Vec3> > temp_mm;
  multimap<int,pair<int,Vec3> > temp_pos_mm;

  // send field id to slaves
  m_comm->broadcast(m_id);


  //receive data from fields
  m_comm->gather(temp_mm);
  m_comm->gather(temp_pos_mm);
  
  // collate receive data from multimap into single map
  for(multimap<int,pair<int,Vec3> >::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    m_save_map.insert(iter->second);
  }
  for(multimap<int,pair<int,Vec3> >::iterator iter=temp_pos_mm.begin();
      iter!=temp_pos_mm.end();
      iter++){
    m_pos_map.insert(iter->second);
  }
}

/*!
  write data out as OpenDX compatible file

  \todo desciption
*/
void VectorParticleFieldMaster::writeAsDX()
{
  //generate filename
  string fn=makeFilename();

  // write header 
  ofstream out_file(fn.c_str());
  out_file << "points = " << m_save_map.size() << endl;
  out_file << "format = ascii" << endl;
  out_file << "dependency = positions, positions" << endl;
  out_file << "interleaving = field" << endl;
  out_file << "field = locations, " << m_field_name << endl;
  out_file << "structure = 3-vector, 3-vector" << endl;
  out_file << "type = float, float  " << endl;
  out_file << "header =  marker \"Start\\n\"" << endl;
  out_file << endl << "end" << endl;
  out_file << "Start" << endl;
  
  // write data
  for(map<int,Vec3>::iterator iter=m_save_map.begin();
      iter!=m_save_map.end();
      iter++){
    out_file << m_pos_map[iter->first] << " " << iter->second << endl;
  }
  out_file.close();
  m_save_map.erase(m_save_map.begin(),m_save_map.end());
  m_pos_map.erase(m_pos_map.begin(),m_pos_map.end());
}

/*!
  write data out as Povray(>=3.0) compatible file
  
  \warning not impl.
*/
void VectorParticleFieldMaster::writeAsPOV()
{
}

/*!
  write data out as SILO file (if supported)
*/
void VectorParticleFieldMaster::writeAsSILO()
{
#if HAVE_SILO
  bool append;
  DBfile* dbfile = openSiloFile(append);

  if (!dbfile) {
    console.Error() << "writeAsSILO not supported";
    return;
  }

  // prepare data
  std::vector<float> x, y, z;
  std::vector<float> vx, vy, vz;
  x.clear(); y.clear(); z.clear();
  vx.clear(); vy.clear(); vz.clear();

  for (map<int,Vec3>::iterator iter=m_save_map.begin();
       iter != m_save_map.end();
       iter++) {
    if (!append) {
      const Vec3& m = m_pos_map[iter->first];
      x.push_back(m.X());
      y.push_back(m.Y());
      z.push_back(m.Z());
    }
    const Vec3& v = iter->second;
    vx.push_back(v.X());
    vy.push_back(v.Y());
    vz.push_back(v.Z());
  }

  // write mesh first if needed
  if (!append) {
    float* coords[] = { &x[0], &y[0], &z[0] };
    DBPutPointmesh(dbfile, "mesh", 3, coords, m_save_map.size(),
                   DB_FLOAT, NULL);
  }

  float* vdata[] = { &vx[0], &vy[0], &vz[0] };
  DBPutPointvar(dbfile, m_field_name.c_str(), "mesh", 3, vdata,
                 m_save_map.size(), DB_FLOAT, NULL);
  // close file 
  DBClose(dbfile);
  m_save_map.erase(m_save_map.begin(),m_save_map.end());
  m_pos_map.erase(m_pos_map.begin(),m_pos_map.end());
#else
  console.Error() << "writeAsSILO not supported\n";
#endif
}

/*!
  sum data and write them out into a single continuous file

  \warning untested
*/
void VectorParticleFieldMaster::writeAsSUM()
{
  Vec3 sum_data=Vec3(0.0,0.0,0.0);
  for(map<int,Vec3>::iterator iter=m_save_map.begin();
      iter!=m_save_map.end();
      iter++){
    sum_data+=iter->second;
  }
  // open file for appending
  ofstream out_file(m_file_name.c_str(),ios::app);
  // write data
  out_file << sum_data << endl;
  // close file
  out_file.close();
  m_save_map.erase(m_save_map.begin(),m_save_map.end());
  m_pos_map.erase(m_pos_map.begin(),m_pos_map.end());
}


/*!
  get max datum and write it out into a single continuous file

  \warning untested
*/
void VectorParticleFieldMaster::writeAsMAX()
{
  Vec3 max_data=(m_save_map.begin())->second;
  for(map<int,Vec3>::iterator iter=m_save_map.begin();
      iter!=m_save_map.end();
      iter++){
    max_data=(iter->second.norm2() > max_data.norm2()) ? iter->second : max_data;
  }
  // open file for appending
  ofstream out_file(m_file_name.c_str(),ios::app);
  // write data
  out_file << max_data << endl;
  // close file
  out_file.close();
  m_save_map.erase(m_save_map.begin(),m_save_map.end());
  m_pos_map.erase(m_pos_map.begin(),m_pos_map.end());
}

/*!
  write data as a raw series of values, one row of values per time step,
  all timesteps into the same file
*/
void VectorParticleFieldMaster::writeAsRAW_SERIES()
{  
  // open file
  ofstream out_file(m_file_name.c_str(),ios::app);
  
  // write data
  for(map<int,Vec3>::iterator iter=m_save_map.begin();
      iter!=m_save_map.end();
      iter++){
    out_file << iter->second << "  ";
  }
  out_file << endl;
  out_file.close();
  m_save_map.erase(m_save_map.begin(),m_save_map.end());
  m_pos_map.erase(m_pos_map.begin(),m_pos_map.end());
}

/*!
  write data as raw position, value pairs, one time step per file
*/
void VectorParticleFieldMaster::writeAsRAW2()
{
  //generate filename
  string fn=makeFilename();

  // open file
  ofstream out_file(fn.c_str());

  // write data
  for(map<int,Vec3>::iterator iter=m_save_map.begin();
      iter!=m_save_map.end();
      iter++){
    out_file << m_pos_map[iter->first] << " " << iter->second << endl;
  }

  // close file
  out_file.close();

  // clean up
  m_save_map.erase(m_save_map.begin(),m_save_map.end());
  m_pos_map.erase(m_pos_map.begin(),m_pos_map.end());
}

/*!
  write data as raw id, position, value triples, one time step per file
*/
void VectorParticleFieldMaster::writeAsRawWithID()
{
  //generate filename
  string fn=makeFilename();

  // open file
  ofstream out_file(fn.c_str());

  // write data
  for(map<int,Vec3>::iterator iter=m_save_map.begin();
      iter!=m_save_map.end();
      iter++){
    out_file << iter->first << " " << m_pos_map[iter->first] << " " << iter->second << endl;
  }

  // close file
  out_file.close();

  // clean up
  m_save_map.erase(m_save_map.begin(),m_save_map.end());
  m_pos_map.erase(m_pos_map.begin(),m_pos_map.end());
}
