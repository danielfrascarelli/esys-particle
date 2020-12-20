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
#include "Fields/WallFieldMaster.h"

// --- STL includes ---
#include <map>

using std::multimap;

// --- IO includes ---
#include <iostream>
#include <fstream>

using std::ofstream;
using std::endl;

/*!
  write the data out as series of raw data, i.e.

  x_0,0 y_0,0 z_0,0 x_1,0 y_1,0 z_1,0...
  x_0,1 y_0,1 z_0,1 x_1,1 y_1,1 z_1,1...
  ...
  where x_i,j is the x-value from wall i at time j
*/
void VectorWallFieldMaster::writeAsRAW_SERIES()
{
  console.XDebug() << "VectorWallFieldMaster::writeAsRAW_SERIES()\n";

  // open file
  ofstream out_file(m_file_name.c_str(),ios::app);

  // set output precision to 10 significant digits, scientific format 
  out_file.setf(ios_base::scientific, ios_base::floatfield);
  out_file.precision(10);

  // write data
  for(map<int,Vec3>::iterator iter=m_data.begin();
      iter!=m_data.end();
      iter++){
    out_file << iter->second << "  ";
  }
  out_file << endl; 
  out_file.close();
  m_data.clear();
}

/*!
  write data out as SILO file (if supported)
*/
void VectorWallFieldMaster::writeAsSILO()
{
#if HAVE_SILO
  bool append;
  DBfile* dbfile = openSiloFile(append);

  if (!dbfile) {
    console.Error() << "writeAsSILO not supported";
    return;
  }

  // prepare data structures
  std::vector<float> x, y, z;
  x.clear(); y.clear(); z.clear();

  for (map<int,Vec3>::iterator iter=m_data.begin();
       iter != m_data.end();
       iter++) {
    const Vec3& v = iter->second;
    x.push_back(v.X());
    y.push_back(v.Y());
    z.push_back(v.Z());
  }

  float* coords[] = { &x[0], &y[0], &z[0] };

  // this is not very nice - the mesh is defined by the Position field
  // which _must_ be written before any other fields
  if (m_field_name == "Position")
    DBPutPointmesh(dbfile, "Position", 3, coords, m_data.size(), DB_FLOAT, NULL);
  else {
    if (!append || DBInqVarExists(dbfile, "Position") == 0)
      console.Error() << "writeAsSILO: You have to write out the Position first\n";
    else
      DBPutPointvar(dbfile, m_field_name.c_str(), "Position", 3, coords,
                     m_data.size(), DB_FLOAT, NULL);
  }
  // close file 
  DBClose(dbfile);
#else
  console.Error() << "writeAsSILO not supported\n";
#endif
  m_data.clear();
}


/*!
  Constructor. Set up the Master and broadcast parameters to the slaves.

  \param comm the communicator
  \param fieldname the name of the field to be saved
  \param filename the name of the file to be saved
  \param walls the names of the walls
  \param savetype the way in which the data are to be saved - only RAW_SERIES supported at the moment
  \param t0 the first timestep to be saved
  \param tend the last timestep to be saved
  \param dt the interval between timesteps to be saving
*/
VectorWallFieldMaster::VectorWallFieldMaster(TML_Comm* comm,
					     const string& fieldname,
					     const string& filename,
					     vector<string> walls,
					     const string& savetype,
					     int t0,
					     int tend,
					     int dt)
  : AFieldMaster(comm,fieldname,filename,savetype,t0,tend,dt)
{
  console.XDebug() << "VectorWallFieldMaster::VectorWallFieldMaster()\n";
  m_comm->broadcast_cont(fieldname);
  console.XDebug() << "bcast fieldname\n";
  m_comm->broadcast(int(walls.size())); 
  for(vector<string>::iterator iter=walls.begin();
      iter!=walls.end();
      iter++){
    m_comm->broadcast_cont(*iter);
  }
  m_comm->broadcast(m_id);
  m_comm->receive(m_sum_flag,1);
  console.XDebug() << "m_sum_flag = " << m_sum_flag << "\n";
  m_comm->barrier();
  console.XDebug() << "end VectorWallFieldMaster::VectorWallFieldMaster()\n";
}

void VectorWallFieldMaster::collect()
{
  console.XDebug() << "VectorWallFieldMaster::collect()\n";
  multimap<int,pair<int,Vec3> > recv_data;

  // send field id to slave
  m_comm->broadcast(m_id);

  // get data from slaves
  m_comm->gather(recv_data);

  if(m_sum_flag==0){ // don't sum -> take data from node 1
    multimap<int,pair<int,Vec3> >::iterator it_lower=recv_data.lower_bound(1);
    multimap<int,pair<int,Vec3> >::iterator it_upper=recv_data.upper_bound(1);
    for(multimap<int,pair<int,Vec3> >::iterator iter=it_lower;
	iter!=it_upper;
	iter++){
      m_data.insert(iter->second);
    }
  } else if(m_sum_flag==1){
    int nslaves=m_comm->size();
    for(int i=1;i<nslaves;i++){
      multimap<int,pair<int,Vec3> >::iterator it_lower=recv_data.lower_bound(i);
      multimap<int,pair<int,Vec3> >::iterator it_upper=recv_data.upper_bound(i);
      for(multimap<int,pair<int,Vec3> >::iterator iter=it_lower;
	  iter!=it_upper;
	  iter++){
	m_data[(iter->second).first]+=(iter->second).second;
      }
    }
  }

  console.XDebug() << "end VectorWallFieldMaster::collect()\n";
}
