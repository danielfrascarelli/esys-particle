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

#include "FieldMaster.h"
#include "field_const.h"

//--- IO includes ---
#include <iostream>
#include <sstream>

using std::cerr;
using std::cout;
using std::endl;
using std::ostringstream;

int AFieldMaster::s_field_count=0;

/*!
  construct AFieldMaster

  \param comm the TML communicator to be used
  \param fieldname the name of the field to be saved
  \param filename the name of the output file
  \param savetype the output file format.  Recognized formats: DX, POV
  \param t0 first timestep to be saved
  \param tend last timestep to be saved

  \param dt save every dt timesteps
*/
AFieldMaster::AFieldMaster(TML_Comm* comm,const string& fieldname,const string& filename,const string& savetype,int t0,int tend,int dt)
{
  cout << "constructing FieldMaster for field " << fieldname << endl;
  m_comm=comm;
  m_field_name=fieldname;
  m_file_name=filename;
  m_t0=t0;
  m_tend=tend;
  m_dt=dt;
  m_save_count=0;
  m_id=s_field_count;
  s_field_count++;
  if(savetype=="DX"){
    m_write_type=WRITE_TYPE_DX;
  } else if(savetype=="POV"){
    m_write_type=WRITE_TYPE_POV;
  } else if(savetype=="SILO"){
    m_write_type=WRITE_TYPE_SILO;
  } else if(savetype=="SUM"){
    m_write_type=WRITE_TYPE_SUM;
  } else if(savetype=="MAX"){
    m_write_type=WRITE_TYPE_MAX;
  } else if(savetype=="RAW_SERIES"){
    m_write_type=WRITE_TYPE_RAW_SERIES;
  } else if(savetype=="RAW2"){
    m_write_type= WRITE_TYPE_RAW2;
  } else if(savetype=="RAW"){
    m_write_type=WRITE_TYPE_RAW;
  } else if(savetype=="RAW_WITH_ID"){
    m_write_type=WRITE_TYPE_RAW_WITH_ID;
  }else if(savetype=="RAW_WITH_POS_ID"){
    m_write_type=WRITE_TYPE_RAW_WITH_POS_ID;
  /****fluid contents: begin****/
  }else if(savetype=="RAW_WITH_POS"){
    m_write_type=WRITE_TYPE_RAW_WITH_POS;
  }else if(savetype=="RAW_WITH_ID_POS"){
    m_write_type=WRITE_TYPE_RAW_WITH_ID_POS;
  }else if(savetype=="RAW_WITH_PARTICLE"){
    m_write_type=WRITE_TYPE_RAW_WITH_PARTICLE;
  }else if(savetype=="VTI"){
    m_write_type=WRITE_TYPE_VTI;
  }else if(savetype=="VTU"){
    m_write_type=WRITE_TYPE_VTU;
  /****fluid contents: end****/
  } else {
    cerr
      << "AFieldMaster: unknown output file format '"
      << savetype
      << "', defaulting to DX" << endl;
  }
}

bool AFieldMaster::needSave(int t)
{
  return (((t-m_t0) % m_dt)==0) && (t>=m_t0) && (t<=m_tend);
}

/*!
  make filename for current save

  \return the filename
*/
string AFieldMaster::makeFilename()
{
  ostringstream fn;
  string suffix;

  switch (m_write_type){
  case WRITE_TYPE_DX : suffix=".dx"; break;
  case WRITE_TYPE_POV : suffix=".pov"; break;
  case WRITE_TYPE_SILO : suffix=".silo"; break;
  case WRITE_TYPE_RAW2 : suffix=".dat"; break;
  case WRITE_TYPE_RAW : suffix=".dat"; break;
  case WRITE_TYPE_RAW_WITH_ID : suffix=".dat"; break;
  case WRITE_TYPE_RAW_WITH_POS_ID : suffix=".dat"; break;
  /****fluid contents: begin****/
  case WRITE_TYPE_RAW_WITH_POS: suffix=".dat"; break; 
  case WRITE_TYPE_RAW_WITH_ID_POS : suffix=".dat"; break;
  case WRITE_TYPE_RAW_WITH_PARTICLE : suffix=".dat"; break; 
  case WRITE_TYPE_VTI : suffix=".vti"; break;
  case WRITE_TYPE_VTU : suffix=".vtu"; break;
  /****fluid contents: end****/
  default : cerr << "AFieldMaster: wrong m_write_type in makeFilename" << endl; // can't happen
  }
  fn << m_file_name << "." << m_save_count*m_dt << suffix;
  m_save_count++;

  return fn.str();
}

#if HAVE_SILO
DBfile* AFieldMaster::openSiloFile(bool& exists)
{
  //generate filename
  string fn = makeFilename();

  DBfile* dbfile = NULL;

  // check if file exists and append data if so
  if (DBInqFile(fn.c_str()) > 0) {
    exists = true;
    // open existing SILO file, try HDF5 format first
    dbfile = DBOpen(fn.c_str(), DB_HDF5, DB_APPEND);
    if (!dbfile)
      dbfile = DBOpen(fn.c_str(), DB_PDB, DB_APPEND);
  } else {
    exists = false;
    // create new SILO file, prefer HDF5 format if available
    dbfile = DBCreate(fn.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_HDF5);
    if (!dbfile)
      dbfile = DBCreate(fn.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
  }
  return dbfile;
}
#endif

/*!
  call the actual write function, depending on m_write_type
*/
void AFieldMaster::write()
{
  switch (m_write_type){
  case WRITE_TYPE_DX : writeAsDX(); break;
  case WRITE_TYPE_POV : writeAsPOV(); break;
  case WRITE_TYPE_SILO : writeAsSILO(); break;
  case WRITE_TYPE_SUM : writeAsSUM(); break;
  case WRITE_TYPE_MAX : writeAsMAX(); break;
  case WRITE_TYPE_RAW_SERIES : writeAsRAW_SERIES(); break;
  case WRITE_TYPE_RAW : writeAsRAW(); break;
  case WRITE_TYPE_RAW2 : writeAsRAW2(); break;
  case WRITE_TYPE_RAW_WITH_ID : writeAsRawWithID(); break;
  case WRITE_TYPE_RAW_WITH_POS_ID : writeAsRawWithPosID(); break;
  /***fluid contents:begin***/
  case WRITE_TYPE_RAW_WITH_POS : writeAsRawWithPos(); break;
  case WRITE_TYPE_RAW_WITH_ID_POS : writeAsRawWithIDPos(); break;
  case WRITE_TYPE_RAW_WITH_PARTICLE : writeAsRawWithParticle(); break;
  case WRITE_TYPE_VTI : writeAsVTI(); break;
  case WRITE_TYPE_VTU : writeAsVTU(); break;
  /***fluid contents:end***/
  default : cerr << "AFieldMaster: wrong m_write_type in write" << endl;
  }
}
