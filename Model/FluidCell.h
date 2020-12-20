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

#ifndef __FLUID_CELL_H
#define __FLUID_CELL_H

// -- project includes --
#include "Foundation/vec3.h"
#include "tml/message/packed_message_interface.h"

/*!
  \class CFluidCell
  \brief base class for fluid cells

  \author Qi Shao
  $Revision$
  $Date$
*/

class CFluidCell
{
 protected:
  double m_Mu;//fluid viscosity
  double m_Phi,m_newPhi,m_effPhi;//porosity
  double m_D, m_K, m_effK;//permeability
  double m_Bf, m_effBf;//fluid bulk modulus
  double m_P, m_disP, m_effP;//pore fluid pressure
  double m_Volume; //total particle volume
  Vec3 m_Size;
  Vec3 m_Vp, m_Vf;//mean particle velocity/fluid velocity
  Vec3 m_Pg;//pressure gradient
  Vec3 m_Pos;//global position of cell centre
  Vec3 m_Index;//global index of cell

  //coefficients for linear equations:
  double c_W;
  double c_E;
  double c_N;
  double c_S;
  double c_D;
  double c_U;
  double c_C;
  double c_B;

 public:
  CFluidCell();
  virtual ~CFluidCell(){};

  inline double getMu() const {return m_Mu;};
  inline double getPhi() const {return m_Phi;};
  inline double getnewPhi() const {return m_newPhi;};
  inline double geteffPhi() const {return m_effPhi;};
  inline double getdPhi() const {return m_newPhi-m_Phi;};
  inline double getD() const {return m_D;};
  inline double getK() const {return m_K;};
  inline double geteffK() const {return m_effK;};
  inline double getBf() const {return m_Bf;};
  inline double geteffBf() const {return m_effBf;};
  inline double getP() const {return m_P;};
  inline double getdisP() const {return m_disP;};
  inline double geteffP() const {return m_effP;};
  inline double getdP() const {return m_disP-m_P;};
  inline double getVolume() const {return m_Volume;};
  inline Vec3 getSize() const {return m_Size;};
  inline Vec3 getVp() const {return m_Vp;};
  inline Vec3 getVf() const {return m_Vf;};
  inline Vec3 getPg() const {return m_Pg;};
  inline Vec3 getPos() const {return m_Pos;};
  inline Vec3 getIndex() const {return m_Index;};
  inline double getAbsVp() const {return m_Vp.norm();};
  inline double getAbsVf() const {return m_Vf.norm();};
  //get coefficients
  inline double getc_W() const {return c_W;};
  inline double getc_E() const {return c_E;};
  inline double getc_N() const {return c_N;};
  inline double getc_S() const {return c_S;};
  inline double getc_D() const {return c_D;};
  inline double getc_U() const {return c_U;};
  inline double getc_C() const {return c_C;};
  inline double getc_B() const {return c_B;};

  void setMu(double Mu){m_Mu=Mu;};
  void setPhi(double Phi){m_Phi=Phi;};
  void setnewPhi(double newPhi){m_newPhi=newPhi;};
  void seteffPhi(double effPhi){m_effPhi=effPhi;};
  void setD(double D){m_D=D;};
  void setK(double K){m_K=K;};
  void seteffK(double effK){m_effK=effK;};
  void setBf(double Bf){m_Bf=Bf;};
  void seteffBf(double effBf){m_effBf=effBf;};
  void setP(double P){m_P=P;};
  void setdisP(double disP){m_disP=disP;};
  void seteffP(double effP){m_effP=effP;};
  void setVolume(double Volume){m_Volume=Volume;};
  void setSize(Vec3 Size){m_Size=Size;};
  void setVp(Vec3 Vp){m_Vp=Vp;};
  void setVf(Vec3 Vf){m_Vf=Vf;};
  void setPg(Vec3 Pg){m_Pg=Pg;};
  void setPos(Vec3 Pos){m_Pos=Pos;};
  void setIndex(Vec3 Index){m_Index=Index;};
  //set coefficients
  void setc_W(double W){c_W=W;};
  void setc_E(double E){c_E=E;};
  void setc_N(double N){c_N=N;};
  void setc_S(double S){c_S=S;};
  void setc_D(double D){c_D=D;};
  void setc_U(double U){c_U=U;};
  void setc_C(double C){c_C=C;};
  void setc_B(double B){c_B=B;};

  CFluidCell& operator=(const CFluidCell&);

  typedef double (CFluidCell::* ScalarFieldFunction)() const;
  typedef Vec3 (CFluidCell::* VectorFieldFunction)() const;

  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);

  friend class TML_PackedMessageInterface;

};


#endif //__FLUID_CELL_H
