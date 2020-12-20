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

#include "Model/ABCDampingIGP.h"

ABCDampingIGP::ABCDampingIGP(const string& type,
			     const string& name,
			     double visc, 
			     double dt,
			     int maxi,
			     const Vec3& vref,
			     const Vec3& pos,
			     const Vec3& normal,
			     double c1)
  :CDampingIGP(type,name,visc,dt,maxi,vref)
{
  m_pos=pos;
  m_normal=normal;
  m_c1=c1;
}

/*!
  Pack the parameters for a ABCDamping into a MPI-buffer

  \param B the buffer
*/
void ABCDampingIGP::packInto(CVarMPIBuffer* B) const
{
  CDampingIGP::packInto(B);
  B->append(m_pos);
  B->append(m_normal);
  B->append(m_c1);
}

/*!
  Extract parameters for a ABCDamping from a MPI-buffer

  \param B the buffer
*/
ABCDampingIGP* extractABCDampingIGP(AMPIBuffer* B)
{
  ABCDampingIGP* res = new ABCDampingIGP;
  res->setName(B->pop_string());
  res->setType(B->pop_string());
  res->setVRef(B->pop_vector());
  res->setVisc(B->pop_double());
  res->setTimeStep(B->pop_double());
  res->setMaxIter(B->pop_int());
  res->setPos(B->pop_vector());
  res->setNormal(B->pop_vector());
  res->setC1(B->pop_double());

  return res;
}
