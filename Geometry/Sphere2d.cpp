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

#include "Sphere2d.h"

double Sphere2D::NearZero=1e-8;

/*!
  Fit circle between 3 circles

  \param iP1 center of 1st circle
  \param iP2 center of 2nd circle
  \param iP3 center of 3rd circle
  \param r1 radius of 1st circle
  \param r2 radius of 2st circle
  \param r3 radius of 3rd circle
  \param M center of fitted circle (output)
  \param r radius of fitted circle (output)
*/
bool Sphere2D::FillIn(const Vec3& P1,const Vec3& P2,const Vec3& P3, double r1, double r2, double r3, Vec3 &M, double &r) 
{
  double B12 = r1*r1-r2*r2+P2*P2-P1*P1 ;
  double B13 = r1*r1-r3*r3+P3*P3-P1*P1 ;
  double V = (cross((P3-P1),(P2-P1))).Z() ;
  if (V==0.0) {
    return false ;
  } 
  Vec3 E = B12*(P3-P1)-B13*(P2-P1) ;
  Vec3 D = (r1-r2)*(P3-P1)-(r1-r3)*(P2-P1) ;
  double a = D.X()*D.X()/(V*V)+ D.Y()*D.Y()/(V*V) -1 ;
  double b = 2.0*((cross(P3,D)).Z())/V + E*D/(V*V) -2.0*r3 ;
  double c = P3*P3+(cross(P3,E)).Z()/V+ E*E/(4*V*V) - r3*r3 ;
  double delta = b*b-4*a*c ;
  if (delta < 0.0) {
    return false ;
  } else if ( delta > 0.0 ) {
    delta = sqrt(delta) ;
  }
  if (a!=0.0) {
    double rs1 = (-b + delta)/(2*a) ;
    double rs2 = (-b - delta)/(2*a) ;
    if ((rs1<=0) && (rs2<=0)) {
      return false ;
    }
    if (rs1<=0) {
      r = rs2 ;
    } else if (rs2 <=0) {
      r = rs1 ;
    } else {
      r = ((rs1)>(rs2) ? (rs2):(rs1)) ;
    }
  } else {
    if (b==0.0) return false ;
    r = -c / b ;
    if (r<=0.0) return false ;
  }
  M.X() = -r*D.Y()/V - E.Y()/(2*V) ;
  M.Y() =  r*D.X()/V + E.X()/(2*V) ;
  return true ;
}

/*!
  Fit circle between 2 circles and a line

  \param iP1 center of 1st circle
  \param iP2 center of 2nd circle
  \param iO origin of line
  \param iD direction of line
  \param r1 radius of 1st circle
  \param r2 radius of 2st circle
  \param M center of fitted circle (output)
  \param r radius of fitted circle (output)
*/
bool Sphere2D::FillInWP(const Vec3& iP1,const Vec3& iP2,const Vec3& iO,const Vec3& iD,double r1,double r2,Vec3& M,double& r) 
{
    Vec3 P2 = iP2-iO ;
    Vec3 P1 = iP1-iO ;
    Vec3 D ;
    D = iD/iD.norm() ;
    Vec3 O = iO-2*(fabs(P1*D)+fabs(P2*D))*D ;
    P2 = iP2-O ;
    P1 = iP1-O ;
    if ((cross(P1,D)).Z()*(cross(P2,D)).Z() < 0.0 ) return false ;
    if ((cross(P1,D)).Z() < 0 ) D *= -1.0 ;	   
    double alpha = 2*P2.X()-2*P1.X()+2*r2*D.Y()-2*r1*D.Y() ;
    double beta = 2*P2.Y()-2*P1.Y()+2*r1*D.X()-2*r2*D.X() ;
    double gamma = P1*P1-P2*P2+r2*r2-r1*r1 ;
    if (gamma == 0.0) {
	return false ;
    }
    if (fabs(beta)<=NearZero) {
        if (alpha == 0.0) return false ;
        double M1x = -gamma/alpha ;
        double a = 1-D.X()*D.X() ;
        double b = 2*D.X()*D.Y()*M1x+2*r1*D.X()-2*P1.Y() ;
        double c = P1*P1-2*P1.X()*M1x+M1x*M1x*(1-D.Y()*D.Y())-2*r1*M1x*D.Y()-r1*r1 ;
        double delta = b*b - 4*a*c ;
        if (delta < 0.0) {
            return false ;
        } else if (delta > 0.0) {
            delta = sqrt(delta) ;
        }
        double rs, My ;
        if (a==0.0) {
            if (b==0.0) return false ;
            My = -c/b ;
            rs = M1x*D.Y()-My*D.X() ;
        } else {
            double M1y = (-b+delta)/(2*a) ;
            double M2y = (-b-delta)/(2*a) ;
            double rs1 = M1x*D.Y()-M1y*D.X() ;
            double rs2 = M1x*D.Y()-M2y*D.X() ;
            if ((rs1 >0) && (rs2>0)) {
                if (rs1<rs2) {
                    rs = rs1 ; My=M1y ;
                } else {
                    rs = rs2 ; My=M2y ;
                }
            } else {
                if (rs1 >0) {
                    rs = rs1 ; My=M1y ;
                } else if (rs2 >0) {
                    rs = rs2 ; My=M2y ;
                } else {
                    return false ;
                }
            }
        }
        M = O+Vec3(M1x,My,0.0) ;        
        r = rs ;
        return true ;
    }
    double A = -alpha/beta ;
    double B = -gamma/beta ;

    double a = A*A+1-(D.Y()-A*D.X())*(D.Y()-A*D.X()) ;
    double b = 2.0*A*B-2*P1.X()-2*P1.Y()*A-2.0*(D.Y()-A*D.X())*(r1-B*D.X()) ;
    double c = P1*P1-2*P1.Y()*B+B*B-(r1-B*D.X())*(r1-B*D.X()) ;

    double delta = b*b - 4*a*c ;
    if (delta < 0.0) {
	return false ;
    } else if (delta > 0.0) {
	delta = sqrt(delta) ;
    }
    double M1x,M2x ;

    if (a!=0) {
	M1x = (-b + delta) / (2*a) ;
	M2x = (-b - delta) / (2*a) ;
    } else {
	if (b==0.0) return false ;
	M1x = -c/b ;
	M2x = M1x ;
    }
    double M1y = A*M1x+B ;
    double M2y = A*M2x+B ;
    double rs1 = M1x*D.Y()-M1y*D.X() ;
    double rs2 = M2x*D.Y()-M2y*D.X() ;
    if ((rs1<=0) && (rs2<=0)) {
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
	M = O+Vec3(M1x,M1y,0.0) ;
    } else {
	M = O+Vec3(M2x,M2y,0.0) ;
    }
    return true;
}

/*!
  Fit circle of fixed radius between circle and line. The equation 
  has always 2 valid solutions, wsol select which solution sould be returned 

  \param iP1 center of circle
  \param iO origin of line
  \param iD direction of line
  \param r1 radius of circle 
  \param M center of fitted circle (output)
  \param r radius of fitted circle (output)
  \param wsol select which solution sould be returned 
*/
bool Sphere2D::FillInWP(const Vec3& iP1,const Vec3& iO,const Vec3& iD,double r1,double r,Vec3& M,int wsol) 
{
  double M1x,M1y ;
  
  Vec3 P1 = iP1-iO ;
  Vec3 D ;
  D = iD/iD.norm() ;
  Vec3 O = iO-2*fabs(P1*D)*D ;
  P1 = iP1-O ;
  if ((cross(P1,D)).Z() < 0 ) D *= -1.0 ;	   
  if (D.X()!=0.0) {
    double a = 1+D.Y()*D.Y()/D.X()/D.X() ;
    double b = -2*D.Y()*r/D.X()/D.X()-2*P1.X()-2*P1.Y()*D.Y()/D.X() ;
    double c = P1*P1+r*r/D.X()/D.X()+2*P1.Y()*r/D.X()-r*r -2*r*r1 -r1*r1;
    if (a==0.0) {
      if (c==0.0) return false ;
      M1x = -c/b ;
      M1y = (r-M1x*D.Y())/D.X() ;
    } else {
      double delta = b*b -4*a*c ;
      double sol = (wsol==1 ? -1.0:1.0) ;
      M1x = (-b + sol*delta)/(2*a) ;
      M1y = (r-M1x*D.Y())/D.X() ;
    }
  } else {
    M1x = r/D.Y() ;
    double a=1 ;
    double b=-2*P1.Y() ;
    double c=P1*P1+M1x*M1x-2*P1.X()*M1x - r*r - 2*r*r1  - r1*r1 ;
    double delta = b*b -4*a*c ;
    double sol = (wsol==1 ? -1.0:1.0) ;
    M1y = (-b + sol*delta)/(2*a) ;
  }
  M = O+Vec3(M1x,M1y,0.0) ;
  return true ;
}
