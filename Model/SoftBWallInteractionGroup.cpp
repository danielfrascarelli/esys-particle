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

#include "Model/SoftBWallInteractionGroup.h"
#include "Foundation/console.h"

//----------------------------------------
//    CSoftBWallIGP member functions 
//----------------------------------------

/*!
  Constructor for bonded wall interaction group with direction dependend elasticity

  \param name the name of the interaction
  \param wallname the name of the wall
  \param kx the spring constant for the elastic interactions in x-direction
  \param ky the spring constant for the elastic interactions in y-direction
  \param kz the spring constant for the elastic interactions in z-direction
  \param tag the tag of the particles to which the wall is bonded (if build via bond and not via distance)
  \param mask the tag mask
  \param scaling toggles scaling of elastic stiffnesses
*/
CSoftBWallIGP::CSoftBWallIGP(const std::string& name,const std::string& wallname,double normalK,double shearK,int tag, int mask, bool scaling)
  : CBWallIGP(name,wallname,normalK,tag,mask)
{
  m_shearK=shearK;
  m_scaling=scaling;
}

void  CSoftBWallIGP::packInto(CVarMPIBuffer* B) const
{
  int sc;
  console.XDebug() << "CSoftBWallIGP::packInto( " << B << " )\n"; 
  CBWallIGP::packInto(B);
  B->append(m_shearK);
  if (m_scaling) {
     sc = 1;
  }
  else {
     sc = 0;
  }
  B->append(sc);
  console.XDebug() << "end CSoftBWallIGP::packInto()\n ";
}

ostream& operator<<(ostream& ost,const CSoftBWallIGP& I)
{
  ost << "CSoftBWallIGP\n";
  ost << "Spring constants : " << I.m_k << " , " << I.m_shearK << endl;
  ost << "Tag             : " << I.m_tag << endl;
  ost << "Scaling         : " << I.m_scaling << endl;
  return ost;
}

CSoftBWallIGP* extractSoftBWallIGP(AMPIBuffer* B)
{
  console.XDebug() << "extractSoftBWallIGP\n";
  string name=B->pop_string();
  double kx=B->pop_double();
  string wallname=B->pop_string();
  int tag=B->pop_int();
  int mask=B->pop_int();
  double ky=B->pop_double();
  int sc = B->pop_int();
  bool scaling;
  if (sc==1) {
     scaling = true;
  }
  else {
     scaling = false;
  }
  CSoftBWallIGP* res=new CSoftBWallIGP(name,wallname,kx,ky,tag,mask,scaling);
  console.XDebug() << "end extractSoftBWallIGP\n";
  return res;
}
