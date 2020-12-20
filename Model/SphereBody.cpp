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

#include "Model/SphereBody.h"
#include "Foundation/console.h"

/*!
  Default constructor. Zeroes all variables. Does _not_ construct
  a useable sphere body!
*/
CSphereBody::CSphereBody()
{
  m_centre=Vec3::ZERO;
  m_radius=0.0;
  m_oldpos=Vec3::ZERO;
  m_force = Vec3::ZERO;
  m_vel   = Vec3::ZERO;
}

/*!
  constructor

  \param c the centre of the sphere body
  \param r the radius of the sphere body
*/
CSphereBody::CSphereBody(const Vec3& c,const double& r)
{
  m_centre=c;
  m_radius=r;
  m_oldpos=m_centre;
  m_force = Vec3::ZERO;
  m_vel   = Vec3::ZERO;
}

/*!
  write restartable checkpoint data to an output stream

  \param ost the output stream
  \param delim
*/
void CSphereBody::writeCheckPoint(ostream& ost,const string &delim) const
{
  ost << m_centre << " " << m_oldpos << " " << m_radius << delim;
}

/*!
  load wall data from a restartable checkpoint

  \param ist the input stream from which the checkpoint is read
*/
void CSphereBody::loadCheckPoint(istream& ist)
{
  ist >> m_centre ;
  ist >> m_oldpos ;
  ist >> m_radius ;
}

ostream& operator<<(ostream& ost,const CSphereBody& w)
{
  ost << "--Wall--" << endl;
  ost << "position : " << w.m_centre << endl;
  ost << "radius   : " << w.m_radius << endl;
  ost << "displ.   : " << w.m_centre-w.m_oldpos << endl;
  ost << flush;
  return ost;
}

