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

#ifndef __MAX_TRIGGER_H
#define __MAX_TRIGGER_H

// --- project includes ---
#include "Foundation/vec3.h"

// --- STL includes ---
#include <map>

using std::map;

struct MaxTrigParams
{
  double trig_on_value;
  double trig_off_value;
  int buff_size;
  int tail_size;
};

class MaxTrigger
{
 private:
  double m_on_level;
  double m_off_level;

 public:
  MaxTrigger(double,double);

  bool Off(const map<int,Vec3>& );
  bool On(const map<int,Vec3>& );
};

#endif // __MAX_TRIGGER_H
