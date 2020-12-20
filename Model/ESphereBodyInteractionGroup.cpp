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

#include "Model/ESphereBodyInteractionGroup.h"
#include "Foundation/console.h"

//----------------------------------------
//    CESphereBodyIGP member functions 
//----------------------------------------

/*!
  Elastic sphere body interaction group constructor

  \param name the name of the sphere body interaction
  \param spherename the name of the sphere body
  \param k the spring constant for the elastic interactions
*/
CESphereBodyIGP::CESphereBodyIGP(
  const  std::string& name,
  const  std::string& spherename,
  double k
)
  : CElasticIGP(name,k)
{
  m_spherename=spherename;
}

void  CESphereBodyIGP::packInto(CVarMPIBuffer* B) const
{
  console.XDebug() << "CESphereBodyIGP::packInto( " << B << " )\n"; 
  CElasticIGP::packInto(B);
  B->append(m_spherename.c_str());
  console.XDebug() << "end CESphereBodyIGP::packInto()\n ";
}

ostream& operator<<(ostream& ost,const CESphereBodyIGP& I)
{
  ost << "CESphereBodyIGP\n";
  ost << "Spring constant : "  << I.m_k << endl;
  return ost;
}

CESphereBodyIGP* extractESphereBodyIGP(AMPIBuffer* B)
{
  string name=B->pop_string();
  double k=B->pop_double();
  string wname=B->pop_string();
  CESphereBodyIGP* res=new CESphereBodyIGP(name,wname,k);
  return res;
}


