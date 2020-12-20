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
#include "VectorTriangleFieldMaster.h"
#include "field_const.h"
#include "console.h"

//--- IO includes ---
#include <iostream>
#include <fstream>
using std::ofstream;
using std::endl;

/*!
  Constructor. Send paramters to slaves

  \param comm the communicator
  \param fieldname the name of the field to be saved
  \param meshname the name of the mesh from which the field is taken
  \param filename the name of the file the field is saved into 
    \param savetype the format in which the data is to be saved
  \param t_0 first timestep to be saved
  \param t_end last timestep to be saved
  \param dt timesteps between saves  
*/
VectorTriangleFieldMaster::VectorTriangleFieldMaster(TML_Comm* comm,
						     const string& fieldname,
						     const string& meshname,
						     const string& filename,
						     const string& savetype,
						     int t_0,int t_end,int dt)
  :AFieldMaster(comm,fieldname,filename,savetype,t_0,t_end,dt)
{
   m_comm->broadcast_cont(fieldname);
   m_comm->broadcast_cont(meshname);
   m_comm->broadcast(m_id);
}

/*!
  Destructor - does nothing
*/
VectorTriangleFieldMaster::~VectorTriangleFieldMaster()
{}
  
/*!
  Call the actual collection function according to the savetype
*/ 
void VectorTriangleFieldMaster::collect()
{
  // send field id to slave
  m_comm->broadcast(m_id);

  if(m_write_type==WRITE_TYPE_DX){
    collectFullDX(); 
  } else {
    collectFull();
  }
}


/*!
  Call the actual writing function according to the savetype
*/ 
void VectorTriangleFieldMaster::write()
{
  switch (m_write_type){
  case WRITE_TYPE_DX: writeAsDX(); break;
  case WRITE_TYPE_RAW:  writeAsRAW(); break;
  default : cerr << "AFieldMaster: wrong m_write_type in write" << endl; 
  }
}

/*!
  write the field in a DX-compatible file
*/
void VectorTriangleFieldMaster::writeAsDX()
{
  console.XDebug() << "VectorTriangleFieldMaster::writeAsDX() - NOT IMPLEMENTED";
}

/*!
  Write the field as raw id,value groups
*/
void VectorTriangleFieldMaster::writeAsRAW()
{
  console.XDebug() << "VectorTriangleFieldMaster::writeAsRAW\n";
  //generate filename
  string fn=makeFilename();

  // open file
  ofstream out_file(fn.c_str());
  // write data
  for(map<int,Vec3>::iterator iter=m_data.begin();
      iter!=m_data.end();
      iter++){
    out_file << iter->first << " " << iter->second << endl;
  }

  // clean up
  out_file.close();
  m_data.clear();
  console.XDebug() << "end VectorTriangleFieldMaster::writeAsRAW\n"; 
}

/*!
  collect the field data from the slaves, just id &
  field value for each triangle
*/
void VectorTriangleFieldMaster::collectFull()
{
  console.XDebug() << "VectorTriangleFieldMaster::collectFull\n";
  multimap<int,pair<int,Vec3> > temp_mm;

  // send type of collect to slaves
  int coll_type=COLL_TYPE_FULL;
  m_comm->broadcast(coll_type);

  // get data from slaves
  m_comm->gather(temp_mm);

  // add data together
  for(multimap<int,pair<int,Vec3> >::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    if(m_data.find((iter->second).first)==m_data.end()){ // id not in m_data -> insert
      m_data.insert(iter->second);
    } else { // id is in m_data -> add
      m_data[(iter->second).first]+=(iter->second).second;
    }
  }
  console.XDebug() << "end VectorTriangleFieldMaster::collectFull\n"; 

}

/*!
  collect the field data from the slaves, including node and
  connectivity info
*/
void VectorTriangleFieldMaster::collectFullDX()
{
  console.XDebug() << "VectorTriangleFieldMaster::collectFullDX() - NOT IMPLEMENTED";

}
