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

#ifndef __FIELDMASTER_H
#define __FIELDMASTER_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if HAVE_SILO
#include <silo.h>
#endif

//--- STL includes ---
#include <string>
using std::string;

// --- project includes ---
#include "Foundation/console.h"


class TML_Comm;

/*!
  \class AFieldMaster
  \brief Abstract base class for master part of field

*/
class AFieldMaster
{
 private:
  static int s_field_count;

 protected:
  TML_Comm *m_comm;
  string m_field_name;
  string m_file_name;
  int m_t0,m_tend,m_dt;
  int m_id;
  int m_save_count;
  int m_write_type; // DX, POV. SILO, SUM, MAX, RAW_SERIES etc

  string makeFilename();
#if HAVE_SILO
  DBfile* openSiloFile(bool& exists);
#endif
  virtual void writeAsDX(){console.Error()<<"writeAsDX NOT IMPLEMENTED\n";};
  virtual void writeAsPOV(){console.Error()<<"writeAsPOV NOT IMPLEMENTED\n";};
  virtual void writeAsSILO(){console.Error()<<"writeAsSILO NOT IMPLEMENTED\n";};
  virtual void writeAsSUM(){console.Error()<<"writeAsSUM NOT IMPLEMENTED\n";};
  virtual void writeAsMAX(){console.Error()<<"writeAsMAX NOT IMPLEMENTED\n";};
  virtual void writeAsRAW_SERIES(){console.Error()<<"writeAsRAW_SERIES NOT IMPLEMENTED\n";};
  virtual void writeAsRAW2(){console.Error()<<"writeAsRAW2 NOT IMPLEMENTED\n";};
  virtual void writeAsRAW(){console.Error()<<"writeAsRAW NOT IMPLEMENTED\n";};
  virtual void writeAsRawWithID(){console.Error()<<"writeAsRawWithID NOT IMPLEMENTED\n";};
  virtual void writeAsRawWithPosID(){console.Error()<<"writeAsRawWithPosID NOT IMPLEMENTED\n";};
  /****fluid contents: begin***/
  virtual void writeAsRawWithPos(){console.Error()<<"writeAsRawWithPosID NOT IMPLEMENTED\n";}; 
  virtual void writeAsRawWithIDPos(){console.Error()<<"writeAsRawWithPosID NOT IMPLEMENTED\n";}; 
  virtual void writeAsRawWithParticle(){console.Error()<<"writeAsRawWithPosID NOT IMPLEMENTED\n";}; 
  virtual void writeAsVTI(){console.Error()<<"writeAsVTI NOT IMPLEMENTED\n";};
  virtual void writeAsVTU(){console.Error()<<"writeAsVTU NOT IMPLEMENTED\n";};
  /****fluid contents: end***/

 public:
  AFieldMaster(TML_Comm*,const string&,const string&,const string&,int,int,int);
  virtual ~AFieldMaster(){};

  virtual bool needSave(int);
  virtual void collect()=0;
  virtual void write();
};

#endif //__FIELDMASTER_H
