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

#include "ViscWallIG.h"
#include "console.h"

//----------------------------------------
//    CVWallIGP member functions 
//----------------------------------------

/*!
  Constructor for viscous wall parameters

  \param name
  \param wallname
  \param k
  \param nu
  \param tag
*/
CVWallIGP::CVWallIGP(const string& name,const string& wallname,double k,double nu,int tag)
  : CEWallIGP(name,wallname,k)
{
  m_tag=tag;
  m_nu=nu;
}

void  CVWallIGP::packInto(CVarMPIBuffer* B) const
{
  console.XDebug() << "CVWallIGP::packInto( " << B << " )\n"; 
  CEWallIGP::packInto(B);
  B->append(m_tag);
  B->append(m_nu);
  console.XDebug() << "end CBWallIGP::packInto()\n ";
}

ostream& operator<<(ostream& ost,const CVWallIGP& I)
{
  ost << "CVWallIGP\n";
  ost << "Spring constant : " << I.m_k << endl;
  ost << "Tag             : " << I.m_tag << endl;
  ost << "Viscosity       : " << I.m_nu << endl;
  return ost;
}

CVWallIGP* extractVWallIGP(AMPIBuffer* B)
{
  console.XDebug() << "extractVWallIGP\n";
  string name=B->pop_string();
  double k=B->pop_double();
  string wname=B->pop_string();
  int tag=B->pop_int();
  double nu=B->pop_double();
  CVWallIGP* res=new CVWallIGP(name,wname,k,nu,tag);
  console.XDebug() << "end extractVWallIGP\n";
  return res;
}
