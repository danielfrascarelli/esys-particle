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


#include "tml/message/packed_message_interface.h"
#include "Model/FluidCell.h" //fluid contents
#include "Foundation/Matrix3.h"

//----
// partial specialisations of the pack/unpack functions
// for simple data types.

//--- int ---
template<>
void TML_PackedMessageInterface::pack<int>(const int& i)
{
  append(i);
}

template<>
void TML_PackedMessageInterface::unpack<int>(int& i)
{
  i=pop_int();
}

//--- double ----
template<>
void TML_PackedMessageInterface::pack<double>(const double& d)
{
  append(d);
}

template<>
void TML_PackedMessageInterface::unpack<double>(double& d)
{
  d=pop_double();
}


/*!
  Pack a Vec3 into a TML packed message

  \param p the Vec3
*/
template<>
void TML_PackedMessageInterface::pack<Vec3>(const Vec3& v)
{
  append(v.X());
  append(v.Y());
  append(v.Z());
}

/*!
  Unpack a Vec3 from a TML packed message

  \param p the Vec3
*/
template<>
void TML_PackedMessageInterface::unpack<Vec3>(Vec3& v)
{
  double db[3];

  pop_doubles(db,3);
  v=Vec3(db[0],db[1],db[2]);
}

/*!
  Pack a Vec3 into a TML packed message

  \param p the Vec3
*/
template<>
void TML_PackedMessageInterface::pack<Matrix3>(const Matrix3& m)
{
    append(m(0,0));
    append(m(0,1));
    append(m(0,2));
    append(m(1,0));
    append(m(1,1));
    append(m(1,2));
    append(m(2,0));
    append(m(2,1));
    append(m(2,2));
}

/*!
  Unpack a Matrix3 from a TML packed message

  \param p the Vec3
*/
template<>
void TML_PackedMessageInterface::unpack<Matrix3>(Matrix3& m)
{
  double db[9];

  pop_doubles(db,9);
  m=Matrix3(db);
}

/***fluid contents: begin****/
/*!
  Pack a CFluidCell into a TML packed message

  \param cell the fluid cell
  */
  template<>
  void TML_PackedMessageInterface::pack<CFluidCell>(const CFluidCell& cell)
  {
    append(cell.m_Mu);
    append(cell.m_Phi);
    append(cell.m_newPhi);
    append(cell.m_effPhi);
    append(cell.m_D);
    append(cell.m_K);
    append(cell.m_effK);
    append(cell.m_Bf);
    append(cell.m_effBf);
    append(cell.m_P);
    append(cell.m_disP);
    append(cell.m_effP);
    append(cell.m_Volume);
    append(cell.m_Size.X());
    append(cell.m_Size.Y());
    append(cell.m_Size.Z());
    append(cell.m_Vp.X());
    append(cell.m_Vp.Y());
    append(cell.m_Vp.Z());
    append(cell.m_Vf.X());
    append(cell.m_Vf.Y());
    append(cell.m_Vf.Z());
    append(cell.m_Pg.X());
    append(cell.m_Pg.Y());
    append(cell.m_Pg.Z());
    append(cell.m_Pos.X());
    append(cell.m_Pos.Y());
    append(cell.m_Pos.Z());
    append(cell.m_Index.X());
    append(cell.m_Index.Y());
    append(cell.m_Index.Z());
    append(cell.c_W);
    append(cell.c_E);
    append(cell.c_N);
    append(cell.c_S);
    append(cell.c_D);
    append(cell.c_U);
    append(cell.c_C);
    append(cell.c_B);
  };


  /*!
    Unpack a CFluidCell from a TML packed message

    \param cell the fluid cell
  */
  template<>
  void TML_PackedMessageInterface::unpack<CFluidCell>(CFluidCell& cell)
  {
    const int numElems = 39;
    double db[numElems] ;

    pop_doubles(db, numElems);
    cell.m_Mu            = db[0];
    cell.m_Phi           = db[1];
    cell.m_newPhi        = db[2];
    cell.m_effPhi        = db[3];
    cell.m_D             = db[4];
    cell.m_K             = db[5];
    cell.m_effK          = db[6];
    cell.m_Bf            = db[7];
    cell.m_effBf         = db[8];
    cell.m_P             = db[9];
    cell.m_disP          = db[10];
    cell.m_effP          = db[11];
    cell.m_Volume        = db[12];
    cell.m_Size          = Vec3(db[13],db[14],db[15]);
    cell.m_Vp            = Vec3(db[16],db[17],db[18]);
    cell.m_Vf            = Vec3(db[19],db[20],db[21]);
    cell.m_Pg            = Vec3(db[22],db[23],db[24]);
    cell.m_Pos           = Vec3(db[25],db[26],db[27]);
    cell.m_Index         = Vec3(db[28],db[29],db[30]);
    cell.c_W             = db[31];
    cell.c_E             = db[32];
    cell.c_N             = db[33];
    cell.c_S             = db[34];
    cell.c_D             = db[35];
    cell.c_U             = db[36];
    cell.c_C             = db[37];
    cell.c_B             = db[38];
  };
  /***fluid contents: end****/
