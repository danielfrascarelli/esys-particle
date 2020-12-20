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

#include "Geometry/Sphere3d.h"

double Sphere3D::NearZero=1e-8;

/*!
  find the sphere that touches 4 spheres 
  
  \param P1 position of 1st Sphere
  \param P2 position of 2nd Sphere
  \param P3 position of 3rd Sphere
  \param P4 position of 4th Sphere
  \param r1 radius of 1st Sphere
  \param r2 radius of 2st Sphere
  \param r3 radius of 3st Sphere
  \param r4 radius of 4st Sphere
  \param M return position of found sphere
  \param r return radius of found sphere
*/
bool Sphere3D::FillIn(const Vec3& P1,const Vec3& P2,const Vec3& P3,const Vec3& P4, double r1, double r2, double r3, double r4, Vec3& M, double& r)
{
  Vec3 A,B,C,D,E,a,b,c,alpha,beta ; 
  double a12,a13,a14,b12,b13,b14 ;

  a.X() = (P2.Y()-P1.Y())*(P4.Y()-P1.Y()) ;
  a.Y() = (P2.Z()-P1.Z())*(P4.Z()-P1.Z()) ;
  a.Z() = (P2.X()-P1.X())*(P4.X()-P1.X()) ;

  b.X() = (P3.Y()-P1.Y())*(P4.Y()-P1.Y()) ;
  b.Y() = (P3.Z()-P1.Z())*(P4.Z()-P1.Z()) ;
  b.Z() = (P3.X()-P1.X())*(P4.X()-P1.X()) ;

  c.X() = (P3.Y()-P1.Y())*(P2.Y()-P1.Y()) ;
  c.Y() = (P3.Z()-P1.Z())*(P2.Z()-P1.Z()) ;
  c.Z() = (P3.X()-P1.X())*(P2.X()-P1.X()) ;

  A.X() = (P3.X()-P1.X())*a.X()-(P2.X()-P1.X())*b.X() ;
  A.Y() = (P3.Y()-P1.Y())*a.Y()-(P2.Y()-P1.Y())*b.Y() ;
  A.Z() = (P3.Z()-P1.Z())*a.Z()-(P2.Z()-P1.Z())*b.Z() ;

  B.X() = (P4.Z()-P1.Z())*c.X()-(P2.Z()-P1.Z())*b.X() ;
  B.Y() = (P4.X()-P1.X())*c.Y()-(P2.X()-P1.X())*b.Y() ;
  B.Z() = (P4.Y()-P1.Y())*c.Z()-(P2.Y()-P1.Y())*b.Z() ;

  C.X() = (P4.X()-P1.X())*c.X()-(P2.X()-P1.X())*b.X() ;
  C.Y() = (P4.Y()-P1.Y())*c.Y()-(P2.Y()-P1.Y())*b.Y() ;
  C.Z() = (P4.Z()-P1.Z())*c.Z()-(P2.Z()-P1.Z())*b.Z() ;

  D.X() = (P3.Z()-P1.Z())*a.X()-(P2.Z()-P1.Z())*b.X() ;
  D.Y() = (P3.X()-P1.X())*a.Y()-(P2.X()-P1.X())*b.Y() ;
  D.Z() = (P3.Y()-P1.Y())*a.Z()-(P2.Y()-P1.Y())*b.Z() ;

  E.X() = A.X()*B.X()-C.X()*D.X() ;
  E.Y() = A.Y()*B.Y()-C.Y()*D.Y() ;
  E.Z() = A.Z()*B.Z()-C.Z()*D.Z() ;

  if ((E.X()==0.0) || (E.Y()==0.0) || (E.Z()==0.0)) {
    return false ;
  }

  a12 = r1-r2 ;
  b12 = r1*r1-r2*r2-P1*P1+P2*P2 ;
  a13 = r1-r3 ;
  b13 = r1*r1-r3*r3-P1*P1+P3*P3 ;
  a14 = r1-r4 ;
  b14 = r1*r1-r4*r4-P1*P1+P4*P4 ;

  alpha.X() = (B.X()*a13*a.X()-B.X()*a12*b.X()-D.X()*a14*c.X()+D.X()*a12*b.X())/E.X() ;
  alpha.Y() = (B.Y()*a13*a.Y()-B.Y()*a12*b.Y()-D.Y()*a14*c.Y()+D.Y()*a12*b.Y())/E.Y() ;
  alpha.Z() = (B.Z()*a13*a.Z()-B.Z()*a12*b.Z()-D.Z()*a14*c.Z()+D.Z()*a12*b.Z())/E.Z() ;
  beta.X() = (B.X()*b13*a.X()-B.X()*b12*b.X()-D.X()*b14*c.X()+D.X()*b12*b.X())/(2.0*E.X()) ;
  beta.Y() = (B.Y()*b13*a.Y()-B.Y()*b12*b.Y()-D.Y()*b14*c.Y()+D.Y()*b12*b.Y())/(2.0*E.Y()) ;
  beta.Z() = (B.Z()*b13*a.Z()-B.Z()*b12*b.Z()-D.Z()*b14*c.Z()+D.Z()*b12*b.Z())/(2.0*E.Z()) ;

  double aa = (alpha*alpha)-1 ;
  double bb = 2.0*(alpha*beta)-2*(alpha*P4)-2*r4 ;
  double cc = (beta*beta)+(P4*P4)-2.0*(beta*P4)-r4*r4 ;
  
  double delta = bb*bb-4*aa*cc ;
  if (delta < 0.0) {
    return false ;
  } else if ( delta > 0.0 ) {
    delta = sqrt(delta) ;
  }
  if (aa!=0.0) {
    double rs1 = (-bb + delta)/(2*aa) ;
    double rs2 = (-bb - delta)/(2*aa) ;
    if ((rs1<=0) && (rs2<=0)) {
      return false ;
    }
    if (rs1<=0) {
      r = rs2 ;
    } else if (rs2 <=0) {
      r = rs1 ;
    } else {
      r = min(rs1,rs2) ;
    }
  } else {
    if (bb==0.0) return false ;
    r = - cc / bb ;
    if (r <=0.0) return false ;
  }
  M = alpha*r+beta ;
  return true ;  
}

/*!
  find the sphere that touch 3 spheres and one wall 
  
  \param P1 position of 1st Sphere
  \param P2 position of 2nd Sphere
  \param P3 position of 3rd Sphere
  \param O origin of the plane
  \param iD normal of the plane
  \param r1 radius of 1st Sphere
  \param r2 radius of 2st Sphere
  \param r3 radius of 3st Sphere
  \param M return position of found sphere
  \param r return radius of found sphere
*/
bool Sphere3D::FillInWP(const Vec3& iP1,const Vec3& iP2,const Vec3& iP3,const Vec3& O,const Vec3& iD, double r1, double r2, double r3, Vec3& M, double& r) 
{
  // cout << "Sphere3D::FillInWP(" << iP1 << ", " << iP2 << ", " << iP3 << ", " << O << ", " << iD << ", " << r1 << ", " << r2 << ", " << r3 << ")" << endl;
    Vec3 P3 = iP3-O ;
    Vec3 P2 = iP2-O ;
    Vec3 P1 = iP1-O ;
    Vec3 W ;
    W = iD/iD.norm() ;
    if ((P1*W) < 0 ) W *= -1.0 ;	   
    if (((P2*W) < 0 ) || ((P3*W) < 0 )){
      // cout << "particles on different sides of wall" << endl;
      return false ;
    }
    Vec3 A12 = 2.0*P2-2.0*P1+(2.0*(r2-r1))*W ;
    Vec3 A13 = 2.0*P3-2.0*P1+(2.0*(r3-r1))*W ;
    double B12 = P1*P1-P2*P2+r2*r2-r1*r1 ;
    double B13 = P1*P1-P3*P3+r3*r3-r1*r1 ;
    double Cx = A12.Z()*A13.Y()-A13.Z()*A12.Y() ;
    double Cy = A12.Y()*A13.X()-A13.Y()*A12.X() ;
    double Cz = A12.Z()*A13.X()-A13.Z()*A12.X() ;
    double Dx = B12*A13.Y()-B13*A12.Y() ;
    double Dy = B12*A13.X()-B13*A12.X() ;
    if (Cy == 0.0) {
      // cout << "Cy == 0.0" << endl;
      return false ;
    }
    Vec3 a = Vec3(Dx/Cy,-Dy/Cy,0.0) ;
    Vec3 b = Vec3(Cx/Cy,-Cz/Cy,1.0) ;
    double aa = b*b-(b*W)*(b*W) ;
    double bb = 2.0*(b*(a-P3))-2.0*(b*W)*(a*W+r3) ;
    double cc = P3*P3-2*(P3*a)+a*a-(a*W+r3)*(a*W+r3) ;
    double delta = bb*bb - 4*aa*cc ;
    if (delta < 0.0) {
      // cout << "delta < 0.0" << endl;
      return false ;
    } else if (delta > 0.0) {
      delta = sqrt(delta) ;
    }
    double M1z,M2z,rs1,rs2 ;
    if (aa!=0) {
      M1z = (-bb + delta) / (2*aa) ;
      M2z = (-bb - delta) / (2*aa) ;
    } else {
      if (bb==0.0) {
        // cout << "bb==0.0" << endl;
        return false ;
      }
      M1z = -cc/bb ;
      M2z = M1z ;
    }
    Vec3 M1,M2 ; 
    M1 = a+b*M1z ;
    M2 = a+b*M2z ;
    rs1 = M1*W ;
    rs2 = M2*W ;
    if ((rs1<=0) && (rs2<=0)) {
      // cout << "rs1,rs2 " << rs1 << " , " << rs2 << endl;
      return false ;
    }
    int sol=1 ;
    if (rs1<=0) {
      r = rs2 ;
      sol=2 ;
    } else if (rs2 <=0) {
      r = rs1 ;
    } else if (rs1<rs2) {
      r = rs1 ;
    } else {
      r = rs2 ;
      sol=2 ;
    }
    if (sol == 1) {
      M = O+M1 ;
    } else {
      M = O+M2 ;
    }
    // cout << "sucess !" << endl;
    return true;
}
