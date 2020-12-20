/////////////////////////////////////////////////////////////
//                                                         //
// Copyright (c) 2003-2014 by The University of Queensland //
// Centre for Geoscience Computing                         //
// http://earth.uq.edu.au/centre-geoscience-computing      //
//                                                         //
// Primary Business: Brisbane, Queensland, Australia       //
// Licensed under the Open Software License version 3.0    //
// http://www.apache.org/licenses/LICENSE-2.0          //
//                                                         //
/////////////////////////////////////////////////////////////

// --- TML includes ---
#include "tml/comm/comm.h"

// --- Project includes ---
#include "FluidInteractionFieldMaster.h"
#include "field_const.h"
#include "console.h"
#include "FieldMaster.h"
#include "Foundation/vec3.h"
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
  Constructor. Setup master and broadcast parameters to slaves

  \param comm the communicator
  \param fieldname the name of the field to be saved
  \param filename the name of the file to be saved into or the base for the generation of the filenames if the saving format requires multiple files
  \param savetype the way to save data, currently supported are DX,SUM
  \param t0 the first timestep to be saved
  \param tend the last timestep to be saved
  \param dt the interval between timesteps to be saving
  \param checked choose between "full" and "checked" fields
*/
ScalarFluidInteractionFieldMaster::ScalarFluidInteractionFieldMaster(TML_Comm* comm,const string& fieldname,const string& filename,const string& savetype,int t0,int tend,int dt)
  :AFieldMaster(comm,fieldname,filename,savetype,t0,tend,dt)
{
  m_comm->broadcast_cont(fieldname);
  m_comm->broadcast(m_id);
}


void ScalarFluidInteractionFieldMaster::collect()
{
  // send field id to slave
  m_comm->broadcast(m_id);

  switch(m_write_type){
  case WRITE_TYPE_SUM: collectSum(); break;
  case WRITE_TYPE_MAX: collectMax(); break;
  case WRITE_TYPE_RAW_WITH_POS: collectFullwithPos(); break;
  case WRITE_TYPE_RAW_WITH_ID: collectFullwithID(); break;
  case WRITE_TYPE_RAW_WITH_ID_POS: collectFullwithIDPos(); break;
  case WRITE_TYPE_RAW_WITH_PARTICLE: collectFullwithParticle(); break;

  default: collectFullwithParticle();
  }
}


/*!
  collect full data set, particle position
*/
void ScalarFluidInteractionFieldMaster::collectFullwithPos()
{
  multimap<int,pair<Vec3,double> > temp_mm;

  // send type of collect to slaves
  int coll_type=COLL_TYPE_FULL_WITH_POS;
  m_comm->broadcast(coll_type);

  // get data from slaves
  m_comm->gather(temp_mm);
  // collate receive data from multimap into single map
  for(multimap<int,pair<Vec3,double> >::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    m_data_with_pos.push_back(iter->second);
  }
}

/*!
  collect full data set, particle id
*/
void ScalarFluidInteractionFieldMaster::collectFullwithID()
{
  multimap<int,pair<int,double> > temp_mm;

  // send type of collect to slaves
  int coll_type=COLL_TYPE_FULL_WITH_ID;
  m_comm->broadcast(coll_type);

  // get data from slaves
  m_comm->gather(temp_mm);
  // collate receive data from multimap into single map
  for(multimap<int,pair<int,double> >::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    m_data_with_id.push_back(iter->second);
  }
}

/*!
  collect full data set, particle id and position
*/
void ScalarFluidInteractionFieldMaster::collectFullwithIDPos()
{
  multimap<int,pair<pair<int,Vec3>,double> > temp_mm;

  // send type of collect to slaves
  int coll_type=COLL_TYPE_FULL_WITH_ID_POS;
  m_comm->broadcast(coll_type);

  // get data from slaves
  m_comm->gather(temp_mm);
  // collate receive data from multimap into single map
  for(multimap<int,pair<pair<int,Vec3>,double> >::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    m_data_with_id_pos.push_back(iter->second);
  }
}

/*!
  collect full data set, particle id, position and radius
*/
void ScalarFluidInteractionFieldMaster::collectFullwithParticle()
{
  multimap<int,pair<esys::lsm::triplet<int,Vec3,double>,double> > temp_mm;

  // send type of collect to slaves
  int coll_type=COLL_TYPE_FULL_WITH_PARTICLE;
  m_comm->broadcast(coll_type);

  // get data from slaves
  m_comm->gather(temp_mm);
  // collate receive data from multimap into single map
  for(multimap<int,pair<esys::lsm::triplet<int,Vec3,double>,double> >::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    m_data_with_particle.push_back(iter->second);
  }
}

/*!
  collect sum of values only
*/
void ScalarFluidInteractionFieldMaster::collectSum()
{
  multimap<int,double> temp_mm;

  // send type of collect to slaves
  int coll_type=COLL_TYPE_SUM;
  m_comm->broadcast(coll_type);

  // get data from slaves
  m_comm->gather(temp_mm);

  // collate receive data from multimap into single map
  for(multimap<int,double>::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    m_sum_vec.push_back(iter->second);
  }
}

/*!
  collect max of values only
*/
void ScalarFluidInteractionFieldMaster::collectMax()
{
  multimap<int,double> temp_mm;

  // send type of collect to slaves
  int coll_type=COLL_TYPE_MAX;
  m_comm->broadcast(coll_type);

  // get data from slaves
  m_comm->gather(temp_mm);

  // collate receive data from multimap into single map
  for(multimap<int,double>::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    m_max_vec.push_back(iter->second);
  }
}


/*!
  write data as pos,value groups
*/
void ScalarFluidInteractionFieldMaster::writeAsRawWithPos()
{
  //generate filename
  string fn=makeFilename();

  // write data
  ofstream out_file(fn.c_str());
  // check if file is sucessfully opened
  if(!out_file){
    console.Error() << "can not open output file " << fn << "\n";
  } else {
    int count=0;
    console.XDebug() << m_data_with_pos.size() << " vectors to be written\n";
    for(vector<pair<Vec3,double> >::iterator iter=m_data_with_pos.begin();
        iter!=m_data_with_pos.end();
        iter++) {
      out_file << iter->first << " " << iter->second << endl;
      count++;
      // write debug message every 10k writes
      if((count % 10000)==0){
        console.XDebug() << count << " vectors written\n";
      }
    }
    console.XDebug() << "finished writing " << count << " vectors \n";
    // close file
    out_file.close();
  }
  //clean up
  m_data_with_pos.clear();
}

/*!
  write data as id,value groups
*/
void ScalarFluidInteractionFieldMaster::writeAsRawWithID()
{
  //generate filename
  string fn=makeFilename();

  // write data
  ofstream out_file(fn.c_str());
  // check if file is sucessfully opened
  if(!out_file){
    console.Error() << "can not open output file " << fn << "\n";
  } else {
    int count=0;
    console.XDebug() << m_data_with_id.size() << " vectors to be written\n";
    for(vector<pair<int,double> >::iterator iter=m_data_with_id.begin();
        iter!=m_data_with_id.end();
        iter++) {
      out_file << iter->first << " " << iter->second << endl;
      count++;
      // write debug message every 10k writes
      if((count % 10000)==0){
        console.XDebug() << count << " vectors written\n";
      }
    }
    console.XDebug() << "finished writing " << count << " vectors \n";
    // close file
    out_file.close();
  }
  //clean up
  m_data_with_id.clear();
}

/*!
  write data as id,pos,value groups
*/
void ScalarFluidInteractionFieldMaster::writeAsRawWithIDPos()
{
  cerr << "ScalarFluidInteractionFieldMaster::writeAsRawWithIDPos not implemented" << endl;
}

/*!
  write data as id,pos,radius,value groups
*/
void ScalarFluidInteractionFieldMaster::writeAsRawWithParticle()
{
  //generate filename
  string fn=makeFilename();

  // write data
  ofstream out_file(fn.c_str());
  // check if file is sucessfully opened
  if(!out_file){
    console.Error() << "can not open output file " << fn << "\n";
  } else {
    int count=0;
    console.XDebug() << m_data_with_particle.size() << " vectors to be written\n";
    for(vector<pair<esys::lsm::triplet<int,Vec3,double>,double> >::iterator iter=m_data_with_particle.begin();
        iter!=m_data_with_particle.end();
        iter++) {
      out_file
        << (iter->first).get<0>()
        << " "
        << (iter->first).get<1>()
        << " "
        << (iter->first).get<2>()
	<< " "
        << iter->second
        << endl;
      count++;
      // write debug message every 10k writes
      if((count % 10000)==0){
        console.XDebug() << count << " vectors written\n";
      }
    }
    console.XDebug() << "finished writing " << count << " vectors \n";
    // close file
    out_file.close();
  }
  //clean up
  m_data_with_particle.clear();
}


/*!
  sum data and write them out into a single continuous file
*/
void ScalarFluidInteractionFieldMaster::writeAsSUM()
{
  // sum data
  double sum_data=0.0;
  for(vector<double>::iterator iter=m_sum_vec.begin();
      iter!=m_sum_vec.end();
      iter++){
    sum_data+=*iter;
  }
  // open file for appending
  ofstream out_file(m_file_name.c_str(),ios::app);

  // set output precision to 10 significant digits, scientific format
  out_file.setf(ios_base::scientific, ios_base::floatfield);
  out_file.precision(10);

  // write data
  out_file << sum_data << endl;

  // close file
  out_file.close();
  //clean up
  m_sum_vec.erase(m_sum_vec.begin(),m_sum_vec.end());

}

/*!
  get the maximum of the data and write it out into a single continuous file
*/
void ScalarFluidInteractionFieldMaster::writeAsMAX()
{
  // sum data
  double max_data=*(m_max_vec.begin());
  for(vector<double>::iterator iter=m_max_vec.begin();
      iter!=m_max_vec.end();
      iter++){
    max_data=(*iter > max_data) ? *iter : max_data;
  }
  // open file for appending
  ofstream out_file(m_file_name.c_str(),ios::app);
  // write data
  out_file << max_data << endl;

  // close file
  out_file.close();
  //clean up
  m_max_vec.erase(m_max_vec.begin(),m_max_vec.end());

}

/*!
  write data as a raw series of values, one row of values per time step,
  all timesteps into the same file
*/
void ScalarFluidInteractionFieldMaster::writeAsRAW_SERIES()
{
  cerr << "ScalarFluidInteractionFieldMaster::writeAsRAW_SERIES not implemented" << endl;
}





// === VECTOR PFM ===
/*!
  Constructor. Setup master and broadcast parameters to slaves

  \param comm the communicator
  \param fieldname the name of the field to be saved
  \param filename the name of the file to be saved into or the base for the generation of the filenames if the saving format requires multiple files
  \param savetype the way to save data, currently supported are DX,SUM
  \param t0 the first timestep to be saved
  \param tend the last timestep to be saved
  \param dt the interval between timesteps to be saving
  \param checked choose between "full" and "checked" fields
*/
VectorFluidInteractionFieldMaster::VectorFluidInteractionFieldMaster(TML_Comm* comm,const string& fieldname,const string& filename,const string& savetype,int t0,int tend,int dt)
  :AFieldMaster(comm,fieldname,filename,savetype,t0,tend,dt)
{
  m_comm->broadcast_cont(fieldname);
  m_comm->broadcast(m_id);
}


void VectorFluidInteractionFieldMaster::collect()
{
  // send field id to slave
  m_comm->broadcast(m_id);

  switch(m_write_type){
  case WRITE_TYPE_SUM: collectSum(); break;
  //case WRITE_TYPE_MAX: collectMax(); break;
  case WRITE_TYPE_RAW_WITH_POS: collectFullwithPos(); break;
  case WRITE_TYPE_RAW_WITH_ID: collectFullwithID(); break;
  case WRITE_TYPE_RAW_WITH_ID_POS: collectFullwithIDPos(); break;
  case WRITE_TYPE_RAW_WITH_PARTICLE: collectFullwithParticle(); break;

  default: collectFullwithParticle();
  }
}


/*!
  collect full data set, particle position
*/
void VectorFluidInteractionFieldMaster::collectFullwithPos()
{
  multimap<int,pair<Vec3,Vec3 > > temp_mm;

  // send type of collect to slaves
  int coll_type=COLL_TYPE_FULL_WITH_POS;
  m_comm->broadcast(coll_type);

  // get data from slaves
  m_comm->gather(temp_mm);
  // collate receive data from multimap into single map
  for(multimap<int,pair<Vec3,Vec3> >::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    m_data_with_pos.push_back(iter->second);
  }
}

/*!
  collect full data set, particle id
*/
void VectorFluidInteractionFieldMaster::collectFullwithID()
{
  multimap<int,pair<int,Vec3> > temp_mm;

  // send type of collect to slaves
  int coll_type=COLL_TYPE_FULL_WITH_ID;
  m_comm->broadcast(coll_type);

  // get data from slaves
  m_comm->gather(temp_mm);
  // collate receive data from multimap into single map
  for(multimap<int,pair<int,Vec3> >::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    m_data_with_id.push_back(iter->second);
  }
}

/*!
  collect full data set, particle id and position
*/
void VectorFluidInteractionFieldMaster::collectFullwithIDPos()
{
  multimap<int,pair<pair<int,Vec3>,Vec3> > temp_mm;

  // send type of collect to slaves
  int coll_type=COLL_TYPE_FULL_WITH_ID_POS;
  m_comm->broadcast(coll_type);

  // get data from slaves
  m_comm->gather(temp_mm);
  // collate receive data from multimap into single map
  for(multimap<int,pair<pair<int,Vec3>,Vec3> >::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    m_data_with_id_pos.push_back(iter->second);
  }
}

/*!
  collect full data set, particle id, position and radius
*/
void VectorFluidInteractionFieldMaster::collectFullwithParticle()
{
  multimap<int,pair<esys::lsm::triplet<int,Vec3,double>,Vec3> > temp_mm;

  // send type of collect to slaves
  int coll_type=COLL_TYPE_FULL_WITH_PARTICLE;
  m_comm->broadcast(coll_type);

  // get data from slaves
  m_comm->gather(temp_mm);
  // collate receive data from multimap into single map
  for(multimap<int,pair<esys::lsm::triplet<int,Vec3,double>,Vec3> >::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    m_data_with_particle.push_back(iter->second);
  }
}

/*!
  collect sum of values only
*/
void VectorFluidInteractionFieldMaster::collectSum()
{
  multimap<int,Vec3> temp_mm;

  // send type of collect to slaves
  int coll_type=COLL_TYPE_SUM;
  m_comm->broadcast(coll_type);

  // get data from slaves
  m_comm->gather(temp_mm);

  // collate receive data from multimap into single map
  for(multimap<int,Vec3>::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    m_sum_vec.push_back(iter->second);
  }
}


/*!
  write data as pos,value groups
*/
void VectorFluidInteractionFieldMaster::writeAsRawWithPos()
{
  //generate filename
  string fn=makeFilename();

  // write data
  ofstream out_file(fn.c_str());
  // check if file is sucessfully opened
  if(!out_file){
    console.Error() << "can not open output file " << fn << "\n";
  } else {
    int count=0;
    console.XDebug() << m_data_with_pos.size() << " vectors to be written\n";
    for(vector<pair<Vec3,Vec3> >::iterator iter=m_data_with_pos.begin();
        iter!=m_data_with_pos.end();
        iter++) {
      out_file << iter->first << " " << iter->second << endl;
      count++;
      // write debug message every 10k writes
      if((count % 10000)==0){
        console.XDebug() << count << " vectors written\n";
      }
    }
    console.XDebug() << "finished writing " << count << " vectors \n";
    // close file
    out_file.close();
  }
  //clean up
  m_data_with_pos.clear();
}

/*!
  write data as id,value groups
*/
void VectorFluidInteractionFieldMaster::writeAsRawWithID()
{
  //generate filename
  string fn=makeFilename();

  // write data
  ofstream out_file(fn.c_str());
  // check if file is sucessfully opened
  if(!out_file){
    console.Error() << "can not open output file " << fn << "\n";
  } else {
    int count=0;
    console.XDebug() << m_data_with_id.size() << " vectors to be written\n";
    for(vector<pair<int,Vec3> >::iterator iter=m_data_with_id.begin();
        iter!=m_data_with_id.end();
        iter++) {
      out_file << iter->first << " " << iter->second << endl;
      count++;
      // write debug message every 10k writes
      if((count % 10000)==0){
        console.XDebug() << count << " vectors written\n";
      }
    }
    console.XDebug() << "finished writing " << count << " vectors \n";
    // close file
    out_file.close();
  }
  //clean up
  m_data_with_id.clear();
}

/*!
  write data as id,pos,value groups
*/
void VectorFluidInteractionFieldMaster::writeAsRawWithIDPos()
{
  cerr << "VectorFluidInteractionFieldMaster::writeAsRawWithIDPos not implemented" << endl;
}

/*!
  write data as id,pos,radius,value groups
*/
void VectorFluidInteractionFieldMaster::writeAsRawWithParticle()
{
  //generate filename
  string fn=makeFilename();

  // write data
  ofstream out_file(fn.c_str());
  // check if file is sucessfully opened
  if(!out_file){
    console.Error() << "can not open output file " << fn << "\n";
  } else {
    int count=0;
    console.XDebug() << m_data_with_particle.size() << " vectors to be written\n";
    for(vector<pair<esys::lsm::triplet<int,Vec3,double>,Vec3> >::iterator iter=m_data_with_particle.begin();
        iter!=m_data_with_particle.end();
        iter++) {
      out_file
        << (iter->first).get<0>()
        << " "
        << (iter->first).get<1>()
        << " "
        << (iter->first).get<2>()
	<< " "
        << iter->second
        << endl;
      count++;
      // write debug message every 10k writes
      if((count % 10000)==0){
        console.XDebug() << count << " vectors written\n";
      }
    }
    console.XDebug() << "finished writing " << count << " vectors \n";
    // close file
    out_file.close();
  }
  //clean up
  m_data_with_particle.clear();
}


/*!
  sum data and write them out into a single continuous file
*/
void VectorFluidInteractionFieldMaster::writeAsSUM()
{
  // sum data
  Vec3 sum_data=Vec3(0.0,0.0,0.0);
  for(vector<Vec3>::iterator iter=m_sum_vec.begin();
      iter!=m_sum_vec.end();
      iter++){
    sum_data+=*iter;
  }
  // open file for appending
  ofstream out_file(m_file_name.c_str(),ios::app);

  // set output precision to 10 significant digits, scientific format
  out_file.setf(ios_base::scientific, ios_base::floatfield);
  out_file.precision(10);

  // write data
  out_file << sum_data << endl;

  // close file
  out_file.close();
  //clean up
  m_sum_vec.erase(m_sum_vec.begin(),m_sum_vec.end());

}

/*!
  write data as a raw series of values, one row of values per time step,
  all timesteps into the same file
*/
void VectorFluidInteractionFieldMaster::writeAsRAW_SERIES()
{
  cerr << "VectorFluidInteractionFieldMaster::writeAsRAW_SERIES not implemented" << endl;
}

