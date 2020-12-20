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

// --- Project includes ---
#include "vec3.h"
#include "Matrix3.h"

/*!
  \class DataParticle
  \brief helper class for a minimal particles just containing the data 
  needed for strain extraction 
*/
class DataParticle
{
 private:
  Vec3 m_pos;
  Vec3 m_initpos;
  double m_rad;
  int m_id;
  // data
  Matrix3 m_tensor_data;
  double m_scalar_data;

 public:
  DataParticle(const Vec3&,const Vec3&,double,int);
  ~DataParticle();

  inline Vec3 getPos() const {return m_pos;};
  inline double getRad() const {return m_rad;};
  inline Vec3 getDisplacement() const{return m_pos-m_initpos;};
  inline int getID() const {return m_id;};
  inline bool isFlagged() const {return false;};
  inline void setFlag(){};
  inline void setTensorData(const Matrix3& M){m_tensor_data=M;};
  inline void setTensorData(int i, int j, double d){m_tensor_data(i,j)=d;};
  inline Matrix3 getTensorData() const {return m_tensor_data;};
  inline double getTensorData(int i,int j){return m_tensor_data(i,j);};
  inline void setScalarData(double d){m_scalar_data=d;};
  inline double getScalarData(){return m_scalar_data;};
};
