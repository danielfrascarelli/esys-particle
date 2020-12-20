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


#include "Model/LocalDampingIGP.h"

//----------------------------------------------
//           LocalDampingIGP functions
//----------------------------------------------

CLocalDampingIGP::CLocalDampingIGP()
  : m_visc(0.0),
    m_dt(0.0)
{
}

CLocalDampingIGP::CLocalDampingIGP(
                         const string& type,
			 const string& name,
			 double        viscosity,
			 double	       dt
)
  : AIGParam(name),
    m_type(type),
    m_visc(viscosity),
    m_dt(dt)
{
}

/*!
  Pack the parameters for a LocalDampingGroup into a MPI-buffer

  \param B the buffer
*/
void CLocalDampingIGP::packInto(CVarMPIBuffer* B) const
{
  AIGParam::packInto(B);
  B->append(m_type.c_str());
  B->append(m_visc);
  B->append(m_dt);
}

/*!
  Extract parameters for a LocalDampingGroup from a MPI-buffer

  \param B the buffer
*/
CLocalDampingIGP* extractLocalDampingIGP(AMPIBuffer* B)
{
  CLocalDampingIGP* res = new CLocalDampingIGP;
  res->setName(B->pop_string());
  res->setType(B->pop_string());
  res->setVisc(B->pop_double());
  res->setTimeStep(B->pop_double());

  return res;
}
