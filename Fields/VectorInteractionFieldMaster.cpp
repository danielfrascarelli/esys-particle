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
#include "VectorInteractionFieldMaster.h"
#include "console.h"
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


/*!
  Constructor. Setup master and broadcast parameters to slaves

  \param comm the communicator
  \param fieldname the name of the field to be saved
  \param igtype the type of interaction group for which the field is saved
  \param igname the name of the interaction group for which the field is saved
  \param filename the name of the file to be saved into or the base for the generation of the filenames if the saving format requires multiple files
  \param savetype the way to save data, currently supported are DX,SUM
  \param t0 the first timestep to be saved
  \param tend the last timestep to be saved
  \param dt the interval between timesteps to be saving
  \param checked choose between "full" and "checked" fields
*/ 
VectorInteractionFieldMaster::VectorInteractionFieldMaster(TML_Comm* comm,const string& fieldname,const string& igtype,const string& igname,const string& filename,const string& savetype,int t0,int tend,int dt,bool checked)
  :AFieldMaster(comm,fieldname,filename,savetype,t0,tend,dt)
{
  m_comm->broadcast_cont(fieldname);
  m_comm->broadcast(m_id);
  m_comm->broadcast_cont(igname);
  m_comm->broadcast_cont(igtype);
  int is_tagged=0;
  m_comm->broadcast(is_tagged);
  int is_checked=checked ? 1 : 0;
  m_comm->broadcast(is_checked);
}

void VectorInteractionFieldMaster::collect()
{
   m_comm->broadcast(m_id);

   switch(m_write_type){
   case WRITE_TYPE_SUM: collectSum(); break;
   case WRITE_TYPE_RAW2: collectFull2();break;
   case WRITE_TYPE_RAW_WITH_ID: collectFullWithID(); break;
   case WRITE_TYPE_RAW_WITH_POS_ID: collectFullWithPosID(); break;
   default: collectFull();
   }
}

/*!
  collect data and <pid1,pid2,pos> info
*/
void VectorInteractionFieldMaster::collectFullWithID()
{
  multimap<int,DataWithID> temp_mm;
  console.XDebug() << "VectorInteractionFieldMaster::collectFullWithID()\n";
  
  // send type of collect to slaves
  int coll_type=COLL_TYPE_FULL_WITH_ID;
  m_comm->broadcast(coll_type);
  // get data from slaves
  m_comm->gather(temp_mm);
  // collate receive data from multimap into single map
  console.XDebug() << temp_mm.size() << " data sets collected\n";
  int count=0;
  for(multimap<int,DataWithID>::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    m_data_with_id.push_back(iter->second);
    count++;
    // write debug message every 10k writes
    if((count % 10000)==0){
      console.XDebug() << count << " data pushed into m_data_with_id\n";
    }  
  }
  console.XDebug() << "finished inserting " << count << " data into m_data_with_id\n";
}

/*!
  collect data and <pid1,pid2,pos1,pos2,ipos> info
*/
void VectorInteractionFieldMaster::collectFullWithPosID()
{
  multimap<int,DataWithPosID> temp_mm;

  // debug output 
  console.XDebug() << "VectorInteractionFieldMaster::collectFullWithPosID()\n";
  
  // send type of collect to slaves
  int coll_type=COLL_TYPE_FULL_WITH_POS_ID;
  m_comm->broadcast(coll_type);
  // get data from slaves
  m_comm->gather(temp_mm);
  // collate receive data from multimap into single map
  console.XDebug() << temp_mm.size() << " data sets collected\n";
  int count=0;
  for(multimap<int,DataWithPosID>::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    m_data_with_pos_id.push_back(iter->second);
    count++;
    // write debug message every 10k writes
    if((count % 10000)==0){
      console.XDebug() << count << " data pushed into m_data_with_id\n";
    }  
  }
  console.XDebug() << "finished inserting " << count << " data into m_data_with_id\n";
}

/*!
  write data out as OpenDX compatible file

  \todo desciption
*/
void VectorInteractionFieldMaster::writeAsDX()
{
  //generate filename
  string fn=makeFilename();

  // write header 
  ofstream out_file(fn.c_str());
  out_file << "points = " << m_data.size() << endl;
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
  for(vector<pair<Vec3,Vec3> >::iterator iter=m_data.begin();
      iter!=m_data.end();
      iter++){
    out_file << iter->first  << " " << iter->second << endl;
  }
  out_file.close(); 
  m_data.erase(m_data.begin(),m_data.end());
}

/*!
  sum data and write them out into a single continuous file
*/
void VectorInteractionFieldMaster::writeAsSUM()
{
  // sum data
  Vec3 sum_data=Vec3(0.0,0.0,0.0);
  for(vector<Vec3>::iterator iter=m_sum_vec.begin();
      iter!=m_sum_vec.end();
      iter++){
    sum_data=sum_data+(*iter);
  }
  // open file for appending
  ofstream out_file(m_file_name.c_str(),ios::app);
  // write data
  out_file << sum_data << endl;

  // close file
  out_file.close();
  //clean up
  m_sum_vec.erase(m_sum_vec.begin(),m_sum_vec.end());

}

/*!
  write data as pos1,pos2,value groups 
*/
void VectorInteractionFieldMaster::writeAsRAW2()
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
    console.XDebug() << m_data2.size() << " vectors to be written\n";
    for(vector<IVecData2>::iterator iter=m_data2.begin();
        iter!=m_data2.end();
        iter++) {
#if 0
      out_file
        << (iter->first).template get<0>()
        << " "
        << (iter->first).template get<1>()
        << " "
        << (iter->first).template get<2>()
        << " "
        << (iter->first).template get<3>()
        << " "
        << (iter->first).template get<4>()
        << " "
        << iter->second
        << endl;
#else
      out_file
        << (iter->first).get<0>()
        << " "
        << (iter->first).get<1>()
        << " "
        << (iter->first).get<2>()
        << " "
        << (iter->first).get<3>()
        << " "
        << (iter->first).get<4>()
        << " "
        << iter->second
        << endl;
#endif
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
  m_data2.clear();
}

/*!
  write data as pid1,pid2,ipos,value groups 
*/
void VectorInteractionFieldMaster::writeAsRawWithID()
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
    for(vector<DataWithID>::iterator iter=m_data_with_id.begin();
        iter!=m_data_with_id.end();
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
  m_data_with_id.clear();
}

/*!
  write data as pid1,pid2,ipos,value groups 
*/
void VectorInteractionFieldMaster::writeAsRawWithPosID()
{
  //generate filename
  string fn=makeFilename();
  
  // debug output
  console.XDebug() << "VectorInteractionFieldMaster::writeAsRawWithPosID() " << fn << "\n"; 

  // write data
  ofstream out_file(fn.c_str());
  // check if file is sucessfully opened
  if(!out_file){
    console.Error() << "can not open output file " << fn << "\n";
  } else {
    int count=0;
    console.XDebug() << m_data_with_id.size() << " vectors to be written\n";
    for(vector<DataWithPosID>::iterator iter=m_data_with_pos_id.begin();
        iter!=m_data_with_pos_id.end();
        iter++) {
      out_file
        << (iter->first).get<0>()
        << " "
        << (iter->first).get<1>()
        << " "
        << (iter->first).get<2>()
	<< " "
	<< (iter->first).get<3>()
        << " "
        << (iter->first).get<4>()
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
  m_data_with_pos_id.clear();
}

/*!
  collect full data set
*/
void VectorInteractionFieldMaster::collectFull()
{
  multimap<int,pair<Vec3,Vec3> > temp_mm;
  
  // send type of collect to slaves
  int coll_type=1;
  m_comm->broadcast(coll_type);

  // get data from slaves
  m_comm->gather(temp_mm);
  // collate receive data from multimap into single map  
  for(multimap<int,pair<Vec3,Vec3> >::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    m_data.push_back(iter->second);
  }
}

/*!
  collect full data set, both particle positions
*/
void VectorInteractionFieldMaster::collectFull2()
{
  multimap<int,IVecData2> temp_mm;

  // send type of collect to slaves
  int coll_type=5;
  m_comm->broadcast(coll_type);
  // get data from slaves
  m_comm->gather(temp_mm);
  // collate receive data from multimap into single map
  console.XDebug() << temp_mm.size() << " data sets collected\n";
  int count=0;
  for(multimap<int,IVecData2>::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    m_data2.push_back(iter->second);
    count++;
    // write debug message every 10k writes
    if((count % 10000)==0){
      console.XDebug() << count << " data pushed into m_data2\n";
    }  
  }
  console.XDebug() << "finished inserting " << count << " data into m_data2\n";
}

/*!
  collect sum of values only 
*/
void VectorInteractionFieldMaster::collectSum()
{
  multimap<int,Vec3> temp_mm;

  // send type of collect to slaves
  m_comm->broadcast(m_write_type);

  // get data from slaves
  m_comm->gather(temp_mm);

  // collate receive data from multimap into single map  
  for(multimap<int,Vec3>::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    m_sum_vec.push_back(iter->second);
  } 
}
