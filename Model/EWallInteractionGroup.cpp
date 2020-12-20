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

#include "Model/EWallInteractionGroup.h"
#include "Foundation/console.h"

//----------------------------------------
//    CEWallIGP member functions 
//----------------------------------------

/*!
  Elastic wall interaction group constructor

  \param name the name of the wall
  \param k the spring constant for the elastic interactions
  \param ipos the initial position of the wall
  \param inorm the initial normal vector of the wall
*/
CEWallIGP::CEWallIGP(
  const  std::string& name,
  const  std::string& wallname,
  double k
)
  : CElasticIGP(name,k)
{
  m_wallname=wallname;
}

void  CEWallIGP::packInto(CVarMPIBuffer* B) const
{
  console.XDebug() << "CEWallIGP::packInto( " << B << " )\n"; 
  CElasticIGP::packInto(B);
  B->append(m_wallname.c_str());
  console.XDebug() << "end CEWallIGP::packInto()\n ";
}

ostream& operator<<(ostream& ost,const CEWallIGP& I)
{
  ost << "CEWallIGP\n";
  ost << "Spring constant : "  << I.m_k << endl;
  return ost;
}

CEWallIGP* extractEWallIGP(AMPIBuffer* B)
{
  string name=B->pop_string();
  double k=B->pop_double();
  string wname=B->pop_string();
  CEWallIGP* res=new CEWallIGP(name,wname,k);
  return res;
}


