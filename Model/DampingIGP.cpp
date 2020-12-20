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


#include "Model/DampingIGP.h"

//----------------------------------------------
//           DampingIGP functions
//----------------------------------------------

CDampingIGP::CDampingIGP()
  : m_vref(Vec3::ZERO),
    m_visc(0.0),
    m_dt(0.0),
    m_max_iter(0)
{
}

CDampingIGP::CDampingIGP(const string& type,
			 const string& name,
			 double        viscosity,
			 double        dt,
			 int           maxIterations,
			 const Vec3&   refVelocity
)
  : AIGParam(name),
    m_type(type),
    m_vref(refVelocity),
    m_visc(viscosity),
    m_dt(dt),
    m_max_iter(maxIterations)
{
}

/*!
  Pack the parameters for a DampingGroup into a MPI-buffer

  \param B the buffer
*/
void CDampingIGP::packInto(CVarMPIBuffer* B) const
{
  AIGParam::packInto(B);
  B->append(m_type.c_str());
  B->append(m_vref);
  B->append(m_visc);
  B->append(m_dt);
  B->append(m_max_iter);
}

/*!
  Extract parameters for a DampingGroup from a MPI-buffer

  \param B the buffer
*/
CDampingIGP* extractDampingIGP(AMPIBuffer* B)
{
  CDampingIGP* res = new CDampingIGP;
  res->setName(B->pop_string());
  res->setType(B->pop_string());
  res->setVRef(B->pop_vector());
  res->setVisc(B->pop_double());
  res->setTimeStep(B->pop_double());
  res->setMaxIter(B->pop_int());

  return res;
}
