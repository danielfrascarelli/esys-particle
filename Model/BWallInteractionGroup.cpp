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

#include "Model/BWallInteractionGroup.h"
#include "Foundation/console.h"

//----------------------------------------
//    CBWallIGP member functions 
//----------------------------------------

/*!
  Bonded wall interaction group constructor

  \param name the name of the interaction
  \param wallname the name of the wall
  \param k the spring constant for the elastic interactions
  \param tag the tag of the particles to which the wall is bonded (if build via bond and not via distance)
  \param mask the particle tag mask
*/
CBWallIGP::CBWallIGP(const std::string& name,const std::string& wallname,double k,int tag, int mask)
  : CEWallIGP(name,wallname,k)
{
  m_tag=tag;
  m_mask=mask;
}

void  CBWallIGP::packInto(CVarMPIBuffer* B) const
{
  console.XDebug() << "CBWallIGP::packInto( " << B << " )\n"; 
  CEWallIGP::packInto(B);
  B->append(m_tag);
  B->append(m_mask);
  console.XDebug() << "end CBWallIGP::packInto()\n ";
}

ostream& operator<<(ostream& ost,const CBWallIGP& I)
{
  ost << "CEWallIGP\n";
  ost << "Name            : " << I.getName() << endl;
  ost << "Wall Name       : " << I.m_wallname << endl;
  ost << "Spring constant : " << I.m_k << endl;
  ost << "Tag             : " << I.m_tag << endl;
  ost << "Mask            : " << I.m_mask << endl;
  return ost;
}

CBWallIGP* extractBWallIGP(AMPIBuffer* B)
{
  console.XDebug() << "extractBWallIGP\n";
  string name=B->pop_string();
  double k=B->pop_double();
  string wallname=B->pop_string();
  int tag=B->pop_int();
  int mask=B->pop_int();
  CBWallIGP* res=new CBWallIGP(name,wallname,k,tag,mask);
  console.XDebug() << "end extractBWallIGP\n";
  return res;
}
