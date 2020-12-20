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

#include "Model/Wall.h"
#include "Foundation/console.h"

/*!
  Default constructor. Zeroes all variables. Does _not_ construct
  a useable wall (normal=(0,0,0)) !
*/
CWall::CWall()
{
  m_origin=Vec3::ZERO;
  m_oldpos=Vec3::ZERO;
  m_normal=Vec3::ZERO;
  m_force = Vec3::ZERO;
  m_vel   = Vec3::ZERO;
}

/*!
  constructor

  \param o the orgin/position of the wall
  \param n the wall normal
*/ 
CWall::CWall(const Vec3& o,const Vec3& n)
{
  m_origin=o;
  m_oldpos=m_origin;
  m_normal=n;
  m_force = Vec3::ZERO;
  m_vel   = Vec3::ZERO;
}

/*!
  generate new vector field slave
  
  \param comm
  \param name
*/
VectorWallFieldSlave<CWall>* CWall::generateVectorFieldSlave(TML_Comm* comm,const string& name)
{
  // get access function
  VectorFieldFunction vf=getVectorFieldFunction(name);

  // create slave
  VectorWallFieldSlave<CWall>* new_fs;
  if(vf!=NULL){
    new_fs=new VectorWallFieldSlave<CWall>(comm,vf);
    new_fs->addWall(this);
  } else {
    new_fs=NULL;   
  }

  return new_fs;
}

/*!
  Get the wall member function which returns a vector field of a given name.

  \param name the name of the field 
*/
CWall::VectorFieldFunction CWall::getVectorFieldFunction(const string& name)
{
  CWall::VectorFieldFunction vf;

  if(name=="Position"){
    vf=&CWall::getPos;
  } else if(name=="Force"){
    vf=&CWall::getForce;
  } else {
    vf=NULL;
    console.Error() << "ERROR - invalid name [ " << name << " ] for wall vector field access function" << "\n"; 
  }

  return vf;
}

/*!
  Get a flag how the field with a given name is to be treated when received by the master, i.e.
  summed over all nodes (Force...) or not (Position...)

  \param fieldname the name of the field
  \returns 1 if the field is to be summed, 0 if not and -1 if the name is invalid
*/
int CWall::getFieldSummationFlag(const string& fieldname)
{
  int res;

  if(fieldname=="Position"){
    res=0;
  } else if(fieldname=="Force"){
    res=1;
  } else {
    res=-1;
    console.Error() << "ERROR - invalid name [ " << fieldname << " ] for wall vector field function" << "\n"; 
  }

  return res;
}


/*!
  write restartable checkpoint data to an output stream

  \param ost the output stream
  \param delim
*/
void CWall::writeCheckPoint(ostream& ost,const string &delim) const
{
  ost << m_origin << " " << m_oldpos << " " << m_normal << delim;
}

/*!
  load wall data from a restartable checkpoint

  \param ist the input stream from which the checkpoint is read
*/
void CWall::loadCheckPoint(istream& ist)
{
  ist >> m_origin ;
  ist >> m_oldpos ;
  ist >> m_normal ;
}

ostream& operator<<(ostream& ost,const CWall& w)
{
  ost << "--Wall--" << endl;
  ost << "position : " << w.m_origin << endl;
  ost << "normal   : " << w.m_normal << endl;
  ost << "displ.   : " << w.m_origin-w.m_oldpos << endl;
  ost << flush;
  return ost;
}

