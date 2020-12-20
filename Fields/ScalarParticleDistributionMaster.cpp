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

#include "tml/comm/comm.h"
#include "ScalarParticleDistributionMaster.h"

/*!
  Constructor without tagging info. Set up the Master and broadcast parameters to the slaves.

  \param comm the communicator
  \param fieldname the name of the field to be saved
  \param filename the name of the file to be saved into or the base for the generation of the filenames if the saving format requires multiple files
  \param savetype the way to save data, currently supported are WINDOW and GLOBAL
  \param t0 the first timestep to be saved
  \param tend the last timestep to be saved
  \param dt_coll the intervall between timesteps for collecting data
  \param dt_save the intervall between timesteps for saving the distribution
  \param x0 minimum value of the field
  \param xmax maximum value of the field
  \param nx number of bins in the histogram
*/
ScalarParticleDistributionMaster::ScalarParticleDistributionMaster(TML_Comm* comm,const string& fieldname,const string& filename,const string& savetype,int t0,int tend,int dt_coll,int dt_save,double x0,double xmax, int nx)
  :ScalarParticleFieldMaster(comm,fieldname,filename,savetype,t0,tend,dt_coll)
{
  m_dt_write=dt_save;
  m_dist=new RealDist(x0,xmax,nx);
  m_is_global=(savetype==string("GLOBAL"));
}

/*!
  Constructor with tagging info. Set up the Master and broadcast parameters to the slaves.

  \param comm the communicator
  \param fieldname the name of the field to be saved
  \param filename the name of the file to be saved into or the base for the generation of the filenames if the saving format requires multiple files
  \param savetype the way to save data, currently supported are WINDOW and GLOBAL
  \param t0 the first timestep to be saved
  \param tend the last timestep to be saved
  \param dt_coll the intervall between timesteps for collecting data
  \param dt_save the intervall between timesteps for saving the distribution
  \param x0 minimum value of the field
  \param xmax maximum value of the field
  \param nx number of bins in the histogram
  \param tag the tag of the particles to be saved
  \param mask the mask to be applied to the tag
*/
ScalarParticleDistributionMaster::ScalarParticleDistributionMaster(TML_Comm* comm,const string& fieldname,const string& filename,const string& savetype,int t0,int tend,int dt_coll,int dt_save,double x0,double xmax, int nx,int tag, int mask)
  :ScalarParticleFieldMaster(comm,fieldname,filename,savetype,t0,tend,dt_coll,tag,mask)
{
  m_dt_write=dt_save;
  m_dist=new RealDist(x0,xmax,nx);
  m_is_global=(savetype==string("GLOBAL"));
}

/*!
  Destructor, deletes the distribution
*/
ScalarParticleDistributionMaster::~ScalarParticleDistributionMaster()
{
  if(m_dist) delete m_dist;
}

/*!
  check if collecting or writing is necessary at current timestep

  \param t the timestep
*/
bool ScalarParticleDistributionMaster::needSave(int t)
{
  bool need_collect;

  need_collect=(((t-m_t0) % m_dt)==0) && (t>=m_t0) && (t<=m_tend);
  m_is_writing_time=(((t-m_t0) % m_dt_write)==0) && (t>=m_t0) && (t<=m_tend);

  return need_collect;
}

/*!
  collect data and add into the distribution
*/
void ScalarParticleDistributionMaster::collect()
{
  // send field id to slave
  m_comm->broadcast(m_id);
 
  // get data from slave
  collectFull();

  // add into distribution
  for(map<int,double>::iterator iter=m_save_map.begin();
      iter!=m_save_map.end();
      iter++){
    m_dist->AddSample(iter->second);
  }
  // clean up
  m_save_map.erase(m_save_map.begin(),m_save_map.end());
  m_pos_map.erase(m_pos_map.begin(),m_pos_map.end());
  m_rad_map.erase(m_rad_map.begin(),m_rad_map.end());
}

/*!
  write data
*/
void ScalarParticleDistributionMaster::write()
{
  if(m_is_writing_time){
    m_dist->Write(m_file_name);
    if(!m_is_global) m_dist->Clear();
  }
}
