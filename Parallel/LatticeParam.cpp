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

#include "Parallel/mpibuf.h"
#include "Parallel/LatticeParam.h"

namespace esys
{
  namespace lsm
  {
    //----------------------------------------------
    //  CLatticeParam functions
    //----------------------------------------------
    CLatticeParam::CLatticeParam(
      const std::string & s,
      double r,
      double a,
      const ProcessDims &dims
    )
    {
      m_particle_type=s;
      m_nrange=r;
      m_alpha=a;
      m_dims = dims;
    }

    void CLatticeParam::packInto(AMPIBuffer *B) const
    {
      B->append(m_nrange);
      B->append(m_alpha);
      B->append(m_particle_type.c_str());
      B->append(static_cast<int>(m_dims[0]));
      B->append(static_cast<int>(m_dims[1]));
      B->append(static_cast<int>(m_dims[2]));
    }

    CLatticeParam CLatticeParam::extractLatticeParam(AMPIBuffer *B)
    {
      double nrange       = B->pop_double();
      double alpha        = B->pop_double();
      std::string particleType = B->pop_string();
      CLatticeParam::ProcessDims dims = CLatticeParam::ProcessDims(3, 0);
      dims[0]             = B->pop_int();
      dims[1]             = B->pop_int();
      dims[2]             = B->pop_int();
    
      return CLatticeParam(particleType, nrange, alpha, dims);
    }

    ostream& operator<<(ostream& ost,const CLatticeParam& CP)
    {
      ost << "CLatticeParam\n";
      ost << "ptype  : " << CP.m_particle_type << endl;
      ost << "nrange : " << CP.m_nrange << endl;
      ost << "alpha  : " << CP.m_alpha << endl;

      return ost;
    }
  }
}
