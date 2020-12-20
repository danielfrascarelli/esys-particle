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

#include "Fields/TriggeredVectorParticleFieldMaster.h"

TriggeredVectorParticleFieldMaster::TriggeredVectorParticleFieldMaster(TML_Comm* comm,
								       const string& fieldname,
								       const string& filename,
								       const string& savetype,
								       int t0,
								       int tend,
								       int dt,
								       const MaxTrigParams& tp) : 
  VectorParticleFieldMaster(comm,fieldname,filename,savetype,t0,tend,dt)
{
  m_save_map_buffer=new RingBuffer<map<int,Vec3> >(tp.buff_size);
  m_pos_map_buffer=new RingBuffer<map<int,Vec3> >(tp.buff_size);
  m_Trigger=new MaxTrigger(tp.trig_on_value,tp.trig_off_value);
  m_tail_size=tp.tail_size;
  m_is_triggered=false;
  m_is_writing_tail=false;
  m_base_file_name=filename;
  m_file_count=0;
  m_ts=0;
}

TriggeredVectorParticleFieldMaster::TriggeredVectorParticleFieldMaster(TML_Comm* comm,
								       const string& fieldname,
								       const string& filename,
								       const string& savetype,
								       int t0,
								       int tend,
								       int dt,
								       int tag,
								       int mask,
								       const MaxTrigParams& tp) : 
  VectorParticleFieldMaster(comm,fieldname,filename,savetype,t0,tend,dt,tag,mask)
{
  m_save_map_buffer=new RingBuffer<map<int,Vec3> >(tp.buff_size);
  m_pos_map_buffer=new RingBuffer<map<int,Vec3> >(tp.buff_size);
  m_Trigger=new MaxTrigger(tp.trig_on_value,tp.trig_off_value);
  m_tail_size=tp.tail_size;
  m_is_triggered=false;
  m_is_writing_tail=false;
  m_base_file_name=filename;
  m_file_count=0;
  m_ts=0;
}

TriggeredVectorParticleFieldMaster::~TriggeredVectorParticleFieldMaster()
{
  delete m_save_map_buffer;
  delete m_pos_map_buffer;
  delete m_Trigger;
}

void TriggeredVectorParticleFieldMaster::write()
{
  m_ts++;
  // check if already triggered
  if(m_is_triggered){ // write
     // check if in tail
    if(m_is_writing_tail){      // if in tail --> write data, decrement tail counter
      // check if still below limit
      if(m_Trigger->On(m_save_map)){ // if above limit -> back to normal writing mode
	m_is_writing_tail=false;
      } else {
	m_tail_counter--;
      }
      if(m_tail_counter==0){ // if end of tail reached --> unset tail flag
	m_is_triggered=false;
	m_is_writing_tail=false;
      }
    } else {  // if not in tail --> check off_trigger
      if(m_Trigger->Off(m_save_map)){ // if off_trigger --> set tail flag, init tail counter, write data
	std::cout << "trigger off at " << m_file_name << m_ts << std::endl;	
	m_is_writing_tail=true;
	m_tail_counter=m_tail_size;
      }
    }
   VectorParticleFieldMaster::write();
  } else {   // if not yet triggered --> check on_trigger
    if(m_Trigger->On(m_save_map)){
       // if on_trigger --> set flag, write out buffer, write data
      IncrementFilename();
      std::cout << "trigger " << m_file_name << " on at  " << m_ts << std::endl;
      m_is_triggered=true;
      m_is_writing_tail=false;
      m_save_map_buffer->insert(m_save_map);
      m_pos_map_buffer->insert(m_pos_map);
      for(int i=0;i<m_save_map_buffer->size();i++){
	m_save_map=(*m_save_map_buffer)[i];
	m_pos_map=(*m_pos_map_buffer)[i];
	VectorParticleFieldMaster::write();
      }
    } else {  // if not --> put data in buffer
      m_save_map_buffer->insert(m_save_map);
      m_pos_map_buffer->insert(m_pos_map);
    }  
  }
  // cleanup maps
  m_save_map.erase(m_save_map.begin(),m_save_map.end());
  m_pos_map.erase(m_pos_map.begin(),m_pos_map.end());
}

void TriggeredVectorParticleFieldMaster::IncrementFilename()
{
  ostringstream fn;

  m_file_count++;
  fn << m_base_file_name << "." << m_file_count;
  m_file_name=fn.str();
}
