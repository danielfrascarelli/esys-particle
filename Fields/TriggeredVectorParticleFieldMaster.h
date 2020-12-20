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

#ifndef __TRIGGERED_VECTOR_PARTICLE_FIELD_MASTER_H
#define __TRIGGERED_VECTOR_PARTICLE_FIELD_MASTER_H

//--- project includes ---
#include "Fields/ParticleFieldMaster.h"
#include "Foundation/RingBuffer.h"
 #include "Fields/MaxTrigger.h"

// --- STL includes ---
#include <map>

using std::map;

class TriggeredVectorParticleFieldMaster : public VectorParticleFieldMaster
{
 private:
  RingBuffer<map<int,Vec3> > *m_save_map_buffer;
  RingBuffer<map<int,Vec3> > *m_pos_map_buffer;
  string m_base_file_name;
  int m_file_count;
  int m_tail_size;
  int m_tail_counter;
  int m_ts;
  bool m_is_triggered;
  bool m_is_writing_tail;

  MaxTrigger* m_Trigger;

  void IncrementFilename();

 public:
  TriggeredVectorParticleFieldMaster(TML_Comm*,const string&,const string&,const string&,int,int,int,const MaxTrigParams&);
  TriggeredVectorParticleFieldMaster(TML_Comm*,const string&,const string&,const string&,int,int,int,int,int,const MaxTrigParams&);
  virtual ~TriggeredVectorParticleFieldMaster();

  virtual void write();
};
#endif // __TRIGGERED_VECTOR_PARTICLE_FIELD_MASTER_H
