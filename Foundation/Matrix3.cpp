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

#include "Foundation/Matrix3.h"
#include "Foundation/Matrix3.hpp"
#include "Foundation/cube_eq.h" // for eigenvalues

#include <iostream>

// --- STL includes ---
#include <set>
#include <cmath>
#include <algorithm>

using std::set;


Matrix3::Matrix3(const double d[9])
{
    m[0][0]=d[0];
    m[0][1]=d[1];
    m[0][2]=d[2];
    m[1][0]=d[3];
    m[1][1]=d[4];
    m[1][2]=d[5];
    m[2][0]=d[6];
    m[2][1]=d[7];
    m[2][2]=d[8];
}

// Solve 3x3 equation system via determinants
// throws exception MatSingularError() in case of a singular matrix
// needs 59 Flops ( 36 mult , 20 add/sub , 3 div )  
Vec3 Matrix3::solve(const Vec3& rhs) const
{
  double det1,detx,dety,detz;

  det1=m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])+m[0][1]*(m[1][2]*m[2][0]-m[1][0]*m[2][2])+m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
  if(det1 == 0){
    throw  MatSingularError();
  }
  detx=rhs.data[0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])+m[0][1]*(m[1][2]*rhs.data[2]-rhs.data[1]*m[2][2])+m[0][2]*(rhs.data[1]*m[2][1]-m[1][1]*rhs.data[2]);
  dety=m[0][0]*(rhs.data[1]*m[2][2]-m[1][2]*rhs.data[2])+rhs.data[0]*(m[1][2]*m[2][0]-m[1][0]*m[2][2])+m[0][2]*(m[1][0]*rhs.data[2]-rhs.data[1]*m[2][0]);
  detz=m[0][0]*(m[1][1]*rhs.data[2]-rhs.data[1]*m[2][1])+m[0][1]*(rhs.data[1]*m[2][0]-m[1][0]*rhs.data[2])+rhs.data[0]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
  double x=detx/det1;
  double y=dety/det1;
  double z=detz/det1;
  return Vec3(x,y,z);
}    

// in situ inversion
// based on a`[ij]=det(A[ij])/det(M)
// where A[ij] is the adjoined submatrix to a[ij]
// 50 Flops (27 mult, 14 add/sub, 9 div)
void Matrix3::invert()
{
  double det1,temp[3][3];

  det1=m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])+m[0][1]*(m[1][2]*m[2][0]-m[1][0]*m[2][2])+m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
  if(det1 == 0){ // det(M)==0 -> singular matrix
    throw  MatSingularError();
  } 
  temp[0][0]=(m[1][1]*m[2][2]-m[1][2]*m[2][1])/det1;
  temp[0][1]=(m[1][2]*m[2][0]-m[1][0]*m[2][2])/det1;
  temp[0][2]=(m[1][0]*m[2][1]-m[1][1]*m[2][0])/det1;
  temp[1][0]=(m[0][2]*m[2][1]-m[0][1]*m[2][2])/det1;
  temp[1][1]=(m[0][0]*m[2][2]-m[0][2]*m[2][0])/det1;
  temp[1][2]=(m[0][1]*m[2][0]-m[0][0]*m[2][1])/det1;
  temp[2][0]=(m[1][2]*m[0][1]-m[1][1]*m[0][2])/det1;
  temp[2][1]=(m[1][0]*m[0][2]-m[1][2]*m[0][0])/det1;
  temp[2][2]=(m[1][1]*m[0][0]-m[1][0]*m[0][1])/det1;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      m[i][j]=temp[i][j];
    }
  }
}

/*!
  Solve the homogeneous system Mx=(0,0,0). It is assumed that rank(M)=2.
  A unit vector is returned;
*/
Vec3 Matrix3::solve_homogeneous() const
{
  Vec3 res;
  
  // check if any colums vector is all zeros
  int idx = 0;
  bool zerocol=false;
  Vec3 cv1=Vec3(m[0][0],m[1][0],m[2][0]);
  Vec3 cv2=Vec3(m[0][1],m[1][1],m[2][1]);
  Vec3 cv3=Vec3(m[0][2],m[1][2],m[2][2]);
  if(cv1.norm()<1e-7){
    idx=0;
    zerocol=true;
  }
  if(cv2.norm()<1e-7){
    idx=1;
    zerocol=true;
  }
  if(cv3.norm()<1e-7){
    idx=2;
    zerocol=true;
  }
  if(zerocol){
    res=Vec3(0.0,0.0,0.0);
    res[idx]=1.0;
  } else {
    // find non-zero element in 1st column
    int i=0;
    while((m[i][0]==0) && (i<3)) i++;
    // the other 2 rows indices
    int j=(i+1)%3;
    int k=(i+2)%3; 

    // copy data
    Matrix3 temp=*this;
    // zero out the 1st column in the other 2 rows
    double f=temp.m[j][0]/temp.m[i][0];
    temp.m[j][0]=0.0; 
    temp.m[j][1]-=f*temp.m[i][1];
    temp.m[j][2]-=f*temp.m[i][2];
    f=temp.m[k][0]/temp.m[i][0];
    temp.m[k][0]=0.0; 
    temp.m[k][1]-=f*temp.m[i][1];
    temp.m[k][2]-=f*temp.m[i][2];
//     cout << "step 1 : " << temp << endl;
    // check if row j is zero (can't have produced a zero column)
    double nj=Vec3(temp.m[j][0],temp.m[j][1],temp.m[j][2]).norm();
    if(nj<1e-7){
      // check if m[k][1] is zero 
      if(fabs(temp.m[k][1])<10e-7){ // diagonalize m[i][0],m[k][2]
    // zero m[i][2]
    f=temp.m[i][2]/temp.m[k][2];
    temp.m[i][2]=0.0;
    temp.m[i][1]-=f*temp.m[k][1];
    //      cout << "step 2a : " << temp << endl;
    // solve for res[0],res[1]
    res[0]=-1.0*temp.m[i][1]/temp.m[i][0];
    res[1]=1.0;
    res[2]=-1.0*temp.m[k][1]/temp.m[k][2];
      } else { // m[k][1]!=0.0 -> diagonalize m[i][0],m[k][1] 
    // zero m[i][1]
    f=temp.m[i][1]/temp.m[k][1];
    temp.m[i][1]=0.0;
    temp.m[i][2]-=f*temp.m[k][2];
    //       cout << "step 2b : " << temp << endl;
    // solve for res[0],res[1]
    res[0]=-1.0*temp.m[i][2]/temp.m[i][0];
    res[1]=1.0;
    res[2]=-1.0*temp.m[k][2]/temp.m[k][1];
      }
    } else {
      // check if m[j][1] is zero 
      if(fabs(temp.m[j][1])<10e-7){ // diagonalize m[i][0],m[j][2]
    // zero m[i][2]
    f=temp.m[i][2]/temp.m[j][2];
    temp.m[i][2]=0.0;
    temp.m[i][1]-=f*temp.m[j][1];
    //      cout << "step 2a : " << temp << endl;
    // solve for res[0],res[1]
    res[0]=-1.0*temp.m[i][1]/temp.m[i][0];
    res[1]=-1.0*temp.m[j][1]/temp.m[j][2];
    res[2]=1.0;
      } else { // m[j][1]!=0.0 -> diagonalize m[i][0],m[j][1] 
    // zero m[i][1]
    f=temp.m[i][1]/temp.m[j][1];
    temp.m[i][1]=0.0;
    temp.m[i][2]-=f*temp.m[j][2];
    //       cout << "step 2b : " << temp << endl;
    // solve for res[0],res[1]
    res[0]=-1.0*temp.m[i][2]/temp.m[i][0];
    res[1]=-1.0*temp.m[j][2]/temp.m[j][1];
    res[2]=1.0;
      }
    }
  }

//   cout << " res: " << res << endl;
  return res.unit();
}

/*!
  calculate eigenvectors & eigenvalues

  \returns v1 eigenvector with largest eigenvalue
  \returns v2 eigenvector with middle eigenvalue
  \returns v3 eigenvector with smalles eigenvalue
  \returns e1 largest eigenvalue
  \returns e2 middle eigenvalue
  \returns e3 smalles eigenvalue
*/
void Matrix3::eigen(Vec3& v1,Vec3& v2,Vec3& v3,double& e1,double& e2,double& e3)
{
  double te1,te2;
  // --- eigenvalues first ---
  // get determinant
  double I1=-1.0*(m[0][0]+m[1][1]+m[2][2]);
  double I2=m[0][0]*m[1][1]-m[0][1]*m[1][0]+
    m[0][0]*m[2][2]-m[0][2]*m[2][0]+
    m[1][1]*m[2][2]-m[1][2]*m[2][1];
  double I3=-1.0*det();
  //std::cout << "I1,I2,I3 " << I1 << " , " << I2 << " , " << I3 << std::endl; 
  if(I3==0){ // simple case, matrix singluar -> 1 eigenvalue 0
    // solve quadratic eqn.
    double p=0.5*I1;
    double q=sqrt(p*p-I2);
    te1=p+q;
    te2=p-q;
    // calc eigenvectors
    Matrix3 m=*this;
    m+=Matrix3::Unit()*-1.0*te1;
    v1=(m.solve(Vec3(0.0,0.0,0.0))).unit();
    m=*this;
    m+=Matrix3::Unit()*-1.0*te2;
    v2=(m.solve(Vec3(0.0,0.0,0.0))).unit();
    v3=cross(v1,v2); 
  } else { // matrix not singular -> solve full cubic equation
    CubicEquation CE(I1,I2,I3);
    set<double> roots=CE.getRealRoots(10e-8);
    if(roots.size()==3){
      set<double>::iterator iter=roots.begin();
      e1=*iter;
      iter++;
      e2=*iter;
      iter++;
      e3=*iter;
      Matrix3 m=*this;
      m+=Matrix3::Unit()*-1.0*e1;
      v1=(m.solve_homogeneous()).unit();
      m=*this;
      m+=Matrix3::Unit()*-1.0*e2;
      v2=(m.solve_homogeneous()).unit();
      m=*this;
      m+=Matrix3::Unit()*-1.0*e3;
      v3=(m.solve_homogeneous()).unit();
    }else{ //can't happen
      // throw some error
    }
  }
}

/*!
    calculate eigenvalues only - not eigenvectors

    \warning expects symmetric matrix
*/
std::array<double,3> Matrix3::eigenvalues() const
{
    std::array<double,3> res;
    double phi;
    const double pi=3.141592653589793;
    
    // check if diagonal 
    double p1=m[0][1]*m[0][1]+m[0][2]*m[0][2]+m[1][2]*m[1][2];
    if (p1==0.0) { // matrix is diagonal -> done, just need sorting
        res[0]=m[0][0];
        res[1]=m[1][1];
        res[2]=m[2][2];
        std::sort(res.begin(),res.end());
        std::swap(res[0],res[2]);
    } else {
        double q=trace()/3.0;
        double h1=m[0][0]-q;
        double h2=m[1][1]-q;
        double h3=m[2][2]-q;
        double p2=h1*h1+h2*h2+h3*h3+2.0*p1;
        double p=sqrt(p2/6.0);
        Matrix3 B=(1.0/p)*(*this-q*Matrix3::Unit());
        double r=B.det()/2.0;
        if(r<=-1.0){
            phi=pi/3.0;
        } else if (r>=1.0) {
            phi = 0.0;
        } else {
            phi = std::acos(r)/3.0;
        }
        
        // get eigenvalues (res[2] -> largest, res[0] -> smallest)
        res[0]=q+2*p*std::cos(phi);
        res[2]=q+2*p*std::cos(phi+(2.0*pi/3.0));
        res[1]=3.0*q-(res[2]+res[0]);
    }
        
    return res;
}

/*!
    extract symmetric part of the matrix, i.e. s_ij=0.5*(m_ij+m_ji)
*/
Matrix3 Matrix3::get_symmetric_part() const
{
    Matrix3 res;
    
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            res.m[i][j]=0.5*(m[i][j]+m[j][i]);
        }
    }
    
    return res;
}

/*!
    calculate the matrix product, i.e. c_ij=sum_k (a_ik * b_kj)
    where 'a' is this matrix

    \param b the other matrix
    \returns the dot product
*/
Matrix3 Matrix3::dot(const Matrix3& b)
{
    Matrix3 res;
    
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            for(int k=0;k<3;k++){
                res(i,j)+=m[i][k]*b.m[k][j];
            }
        }
    }
    
    return res;
}