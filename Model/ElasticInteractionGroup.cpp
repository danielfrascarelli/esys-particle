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

#include "ElasticInteractionGroup.h"
#include "console.h"

#include <memory>

//----------------------------------------
//    CElasticIGP member functions 
//----------------------------------------

void CElasticIGP::packInto(CVarMPIBuffer* B) const
{
  AIGParam::packInto(B);
  B->append(m_k);
}

ostream& operator<<(ostream& ost,const CElasticIGP& B)
{
  ost << "CElasticIGP\n";
  ost << "Spring Constant : " << B.m_k << endl;
  return ost;
}
