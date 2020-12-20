/////////////////////////////////////////////////////////////
//                                                         //
// Copyright (c) 2003-2014 by The University of Queensland //
// Centre for Geoscience Computing                         //
// http://earth.uq.edu.au/centre-geoscience-computing      //
//                                                         //
// Primary Business: Brisbane, Queensland, Australia       //
// Licensed under the Open Software License version 3.0    //
// http://www.apache.org/licenses/LICENSE-2.0          //
//                                                         //
/////////////////////////////////////////////////////////////


// -- project includes --
#include "Model/FluidCell.h"
#include "tml/message/packed_message_interface.h"


CFluidCell::CFluidCell() 
{
  m_Mu=0;m_Phi=0;m_newPhi=0;m_effPhi=0;m_D=0;m_K=0;m_effK=0;m_Bf=0;m_effBf=0;m_P=0;m_disP=0;m_effP=0;m_Volume=0;m_Size=Vec3(0,0,0);
  m_Vp=Vec3(0,0,0);m_Vf=Vec3(0,0,0);m_Pg=Vec3(0,0,0);m_Pos=Vec3(0,0,0);
  c_W=0;c_E=0;c_N=0;c_S=0;c_D=0;c_U=0;c_C=0;c_B=0;
}


CFluidCell& CFluidCell::operator=(const CFluidCell& rhs)
{
  m_Mu=rhs.m_Mu;
  m_Phi=rhs.m_Phi;
  m_newPhi=rhs.m_newPhi;
  m_effPhi=rhs.m_effPhi;
  m_D=rhs.m_D;
  m_K=rhs.m_K;
  m_effK=rhs.m_effK;
  m_Bf=rhs.m_Bf;
  m_effBf=rhs.m_effBf;
  m_P=rhs.m_P;
  m_disP=rhs.m_disP;
  m_effP=rhs.m_effP;
  m_Volume=rhs.m_Volume;
  m_Size=rhs.m_Size;
  m_Vp=rhs.m_Vp;
  m_Vf=rhs.m_Vf;
  m_Pg=rhs.m_Pg;
  m_Pos=rhs.m_Pos;
  m_Index=rhs.m_Index;
  c_W=rhs.c_W;
  c_E=rhs.c_E;
  c_N=rhs.c_N;
  c_S=rhs.c_S;
  c_D=rhs.c_D;
  c_U=rhs.c_U;
  c_C=rhs.c_C;
  c_B=rhs.c_B;
  return *this;
}


/*!
  Get the fluidcell member function which returns a scalar field of a given name.

  \param name the name of the field
*/
CFluidCell::ScalarFieldFunction CFluidCell::getScalarFieldFunction(const string& name)
{
  CFluidCell::ScalarFieldFunction sf;

  if(name=="Mu"){
    sf=&CFluidCell::getMu;
  } else if(name=="Phi"){
    sf=&CFluidCell::geteffPhi;
  } else if(name=="dPhi"){
    sf=&CFluidCell::getdPhi;
  } else if(name=="D"){
    sf=&CFluidCell::getD;
  } else if(name=="K"){
    sf=&CFluidCell::geteffK;
  } else if(name=="Bf"){
    sf=&CFluidCell::geteffBf;
  } else if(name=="P"){
    sf=&CFluidCell::getP;
  } else if(name=="disP"){
    sf=&CFluidCell::getdisP;
  } else if(name=="effP"){
    sf=&CFluidCell::geteffP;
  } else if(name=="dP"){
    sf=&CFluidCell::getdP;
  } else if(name=="Volume"){
    sf=&CFluidCell::getVolume;
  } else if(name=="Vp_abs"){
    sf=&CFluidCell::getAbsVp;
  } else if(name=="Vf_abs"){
    sf=&CFluidCell::getAbsVf;
  } else if(name=="c_W"){
    sf=&CFluidCell::getc_W;
  } else if(name=="c_E"){
    sf=&CFluidCell::getc_E;
  } else if(name=="c_N"){
    sf=&CFluidCell::getc_N;
  } else if(name=="c_S"){
    sf=&CFluidCell::getc_S;
  } else if(name=="c_D"){
    sf=&CFluidCell::getc_D;
  } else if(name=="c_U"){
    sf=&CFluidCell::getc_U;
  } else if(name=="c_C"){
    sf=&CFluidCell::getc_C;
  } else if(name=="c_B"){
    sf=&CFluidCell::getc_B;
  } else {
    sf=NULL;
    cerr << "ERROR - invalid name for fluid cell scalar  access function" << endl;
  }

  return sf;
}


/*!
  Get the fluid cell member function which returns a vector field of a given name.

  \param name the name of the field
*/

CFluidCell::VectorFieldFunction CFluidCell::getVectorFieldFunction(const string& name)
  {
  CFluidCell::VectorFieldFunction sf;

  if(name=="Vp"){
    sf=&CFluidCell::getVp;
  } else if(name=="Vf"){
    sf=&CFluidCell::getVf;
  } else if(name=="Pos"){
    sf=&CFluidCell::getPos;
  } else if(name=="Index"){
    sf=&CFluidCell::getIndex;
  } else {
    sf=NULL;
    cerr << "ERROR - invalid name for fluid cell vector access function" << endl;
  }

  return sf;
}

