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

#include <mpi.h>
#include "Model/RotThermBondedInteraction.h"
#include "Model/BondedInteractionCpData.h"
#include "Foundation/console.h"

// --- TML includes ---
#include "tml/message/packed_message_interface.h"

#include <stdexcept>
/*
double calc_angle( double s_in, double c_os)
{
  double angle;
  if (s_in >0.0)
  {
    if ( c_os>=1.0 ||c_os<=-1.0)
    {
      angle = 0.0;
    }
    else {
      angle  =   acos(c_os);
    }
  }
  else if (s_in==0.0)
  {
    angle = 0.0;      
  }
  else
  {
    if ( c_os>=1.0 ||c_os<=-1.0)
    {
      angle = 0.0;
    }
    else
    {
      angle  = - acos(c_os);
    }
  }
  return   angle ;
}
*/

CRotThermBondedIGP::CRotThermBondedIGP()
 : kr(0.0),
   ks(0.0),
   kt(0.0),
   kb(0.0),
   max_nForce(0.0),
   max_shForce(0.0),
   max_tMoment(0.0),
   max_bMoment(0.0),
   diffusivity(0.0),
   tag(0)
{
}

CRotThermBondedIGP::CRotThermBondedIGP(
    const std::string &name,
    double normalK,
    double shearK,
    double torsionK,
    double bendingK,
    double maxNormalForce,
    double maxShearForce,
    double maxTorsionMoment,
    double maxBendingMoment,
    double thermalDiffusivity,
    int bondTag
)
 : AIGParam(name),
   kr(normalK),
   ks(shearK),
   kt(torsionK),
   kb(bendingK),
   max_nForce(maxNormalForce),
   max_shForce(maxShearForce),
   max_tMoment(maxTorsionMoment),
   max_bMoment(maxBendingMoment),
   diffusivity(thermalDiffusivity),
   tag(bondTag)
{
}

CRotThermBondedInteraction::CRotThermBondedInteraction():ARotThermPairInteraction()
{
  m_kr = 0.0 ;
  m_ks = 0.0 ;
  m_kb = 0.0 ;
  m_kt = 0.0 ;

  m_max_nForce  = 0.0 ;
  m_max_shForce = 0.0 ;
  m_max_tMoment = 0.0 ;
  m_max_bMoment = 0.0 ;

  m_nForce  = 0.0 ;
  m_shForce = 0.0;
  m_tMoment = 0.0;
  m_bMoment = 0.0;
  m_min_r =0.0  ;
  m_moment  = Vec3(0.0,0.0,0.0);
  m_force   = Vec3(0.0,0.0,0.0);

  m_diffusivity =0.0;

  m_tag = 0;
}

CRotThermBondedInteraction::CRotThermBondedInteraction(CRotThermParticle* p1,CRotThermParticle* p2,const CRotThermBondedIGP& param):ARotThermPairInteraction(p1,p2)
{
  m_nForce  = 0.0;
  m_shForce = 0.0;
  m_tMoment = 0.0;
  m_bMoment = 0.0;

//12/02/2005  added

  double ran_ratio ;
  double ran_ratioH ;
  double km_coef ;

  if (m_p1->getRad()<= m_p2->getRad() )    m_min_r = m_p1->getRad() ;
  else                                     m_min_r = m_p2->getRad() ;


  if(m_p1->getDo2dCalculations()) { // 2D
    ran_ratio = 2.0*m_min_r/(m_p1->getRad()+m_p2->getRad());
    km_coef = m_min_r*m_min_r;
    ran_ratioH = 2.0*m_min_r*(m_p1->getRad()+m_p2->getRad());
  }else{ //  3D
    ran_ratio = 2.0*m_min_r*m_min_r/(m_p1->getRad()+m_p2->getRad());
   km_coef = 1.0;  //tmp ????
   ran_ratioH = 2.0*m_min_r*m_min_r*(m_p1->getRad()+m_p2->getRad());
  }

 
  m_kr =  ran_ratio*param.kr  ;
  m_ks =  ran_ratio*param.ks  ;
  m_kb =  km_coef*ran_ratio*param.kb  ;
  m_kt =  ran_ratio*param.kt  ;  // temp ??


//  cout << " m_kr= " <<m_kr<<" m_ks =  " <<m_ks <<" m_kb=  " <<m_kb <<" m_kt=  " <<m_kt << endl;





  //  double rand_num = (double) (rand())/(double)(RAND_MAX) ;

  //rand_num = 0.65+ 0.35*rand_num ;
  double  rand_num = 1.0;


  if(m_p1->getDo2dCalculations()) {
  m_max_nForce =  m_min_r*param.max_nForce*rand_num ;
  m_max_shForce = m_min_r*param.max_shForce*rand_num ;
  m_max_tMoment = m_min_r*param.max_tMoment*rand_num ;
  m_max_bMoment = km_coef*m_min_r*param.max_bMoment*rand_num ;

//  cout << "  m_max_nForce= " <<m_max_nForce<< "  m_max_shForce= " <<m_max_shForce <<"  m_max_tMoment= " <<m_max_tMoment << "  m_max_bMoment=  " <<m_max_bMoment << endl;

  }else{

  m_max_nForce =  m_min_r*m_min_r*param.max_nForce*rand_num ;
  m_max_shForce = m_min_r*m_min_r*param.max_shForce*rand_num ;
  m_max_tMoment = m_min_r*m_min_r*param.max_tMoment*rand_num ;
  m_max_bMoment = m_min_r*m_min_r*param.max_bMoment*rand_num ;
 //cout << "  m_max_nForce= " <<m_max_nForce<< "  m_max_shForce= " <<m_max_shForce <<"  m_max_tMoment= " <<m_max_tMoment << "  m_max_bMoment=  " <<m_max_bMoment << endl;

  }

  m_diffusivity = ran_ratioH*param.diffusivity ; // ratio ???


 //24/Mar, 2006    strong boundary
/*
  if (m_p1->getTag()*m_p2->getTag()!=0 ) {

   m_kr *= 100.0 ;
   m_ks *= 100.0 ;
   m_kb *= 100.0 ;
   m_kt *= 100.0 ;
   m_max_nForce *= 100.0 ;
   m_max_shForce*= 100.0 ;
   m_max_tMoment*= 100.0 ;
   m_max_bMoment*= 100.0 ;
  }
*/



// cout <<m_p1->getID() <<" -----   " << m_p2->getID() << endl;
// cout << "ran_ratio=  " <<ran_ratio << endl;
// cout << "m_kr  " << m_kr <<"   " <<  m_ks << "  " <<m_kb << "  "  << m_kt << endl;

// cout << "m_max_nForce= " <<m_max_nForce <<"   "  << m_max_shForce<<"   " <<m_max_tMoment <<"   " <<m_max_bMoment << endl;


/*

  m_kr =  param.kr  ;
  m_ks =  param.ks  ;
  m_kb =  param.kb  ;
  m_kt =  param.kt  ;

  m_max_nForce =  param.max_nForce ;
  m_max_shForce = param.max_shForce ;
  m_max_tMoment = param.max_tMoment ;
  m_max_bMoment = param.max_bMoment ;



*/


  m_force  = Vec3(0.0,0.0,0.0);
  m_moment = Vec3(0.0,0.0,0.0) ;

  const Vec3 D = p1->getPos()-p2->getPos();
  m_dist=sqrt(D*D);
  
  m_tag = param.tag;
}

int CRotThermBondedInteraction::getTag() const
{
  return m_tag;
}

void CRotThermBondedInteraction::setTag(int tag)
{
  m_tag = tag;
}



CRotThermBondedInteraction::~CRotThermBondedInteraction()
{
}

/*
bool CRotBondedInteraction::broken()
{
  const double crit_nf = (m_nForce/m_max_nForce > 0.0) ? m_nForce/m_max_nForce :0.0;
  const double criterion =
    crit_nf  +
    m_shForce/m_max_shForce +
    m_tMoment/m_max_tMoment +
    m_bMoment/m_max_bMoment;

  // if(criterion > 1.0 ) {
//     std::stringstream msg;
//     msg << "bond broken" << "\n" ;
//     msg << "ids : " << m_p1->getID() <<  " " << m_p2->getID() << "\n" ;
//     msg << "positions : " << m_p1->getPos() <<  m_p2->getPos() << "\n";
//     msg << "dist : " << m_dist << "\n";
//     msg << "m_nForce/m_max_nForce = "   << m_nForce/m_max_nForce   << "\n";
//     msg << "m_shForce/m_max_shForce = " << m_shForce/m_max_shForce << "\n";
//     msg << "m_tMoment/m_max_tMoment = " << m_tMoment/m_max_tMoment << "\n";
//     msg << "m_bMoment/m_max_bMoment = " << m_bMoment/m_max_bMoment << "\n";
//     console.Debug() << msg.str();
//   }
  return(criterion > 1.0);
}

*/


bool CRotThermBondedInteraction::broken()
{

  bool res;

// const double crit_nf = (m_nForce/m_max_nForce > 0.0) ? m_nForce/m_max_nForce :0.0;
//  const double criterion =
//    crit_nf  +
//    m_shForce/m_max_shForce +
//    m_tMoment/m_max_tMoment +
//    m_bMoment/m_max_bMoment;



  double criterion = m_nForce/m_max_nForce +
                     m_shForce*m_shForce/(m_max_shForce*m_max_shForce) +
                     m_tMoment/m_max_tMoment +
                     m_bMoment/m_max_bMoment    ;
//         cout << "wyc crit= " << criterion<< endl;
//         cout << "m_nForce= " <<m_nForce <<"  " << m_max_nForce <<endl;
//         cout << "m_shForce= " <<m_shForce <<"  " << m_max_shForce <<endl;
//         cout << "m_tMoment= " <<m_tMoment <<"  " << m_max_tMoment <<endl;
//         cout << "m_bMoment= " <<m_bMoment <<"  " << m_max_bMoment <<endl;


 //  int shear_tensile = 0  ;

  if(criterion > 1.0 ){
    console.Debug() << "bond broken" << "\n" ;
    console.Debug() << "ids : " << m_p1->getID() <<  " " << m_p2->getID() << "\n" ;
    console.Debug() << "positions : " << m_p1->getPos() <<  m_p2->getPos() << "\n";
    console.Debug() << "dist : " << m_dist << "\n" ;

    if(m_p1!=NULL) m_p1->setFlag();
    if(m_p2!=NULL) m_p2->setFlag(); 
    res=true;

    //console.Debug() << "break : " << m_break << "\n" ;

//zjr    
//      cout << "wyc bond broken " << m_p1->getID()<< " ---   "<<  m_p2->getID() <<endl;

//      cout << "wyc  Rad " << m_p1->getRad()<< " ---   "<<  m_p2->getRad() <<endl;
//      cout << " jixiang  "<<m_nForce/m_max_nForce<< endl;
//      cout << "qie xiang  " << m_shForce/m_max_shForce << endl;
//      cout << " wan qu  "  << m_bMoment/m_max_bMoment << endl;
//      cout << "niu zuan  " <<  m_tMoment/m_max_tMoment << endl;
//      cout << " huai  "<<  (m_p1->getOldPos() + m_p2->getOldPos())*0.5<< endl; 
//zjr 
//       cout << "wyc bond broken " << m_p1->getPos()<< " ---   "<<  m_p2->getPos() <<endl;
//zjr  
//      cout << "wyc bond broken " << m_p1->getQuat()<< " ---   "<<  m_p2->getQuat() <<endl;

   /*
    if (m_nForce<0.0) {

    cout << "m_nForce= " <<m_nForce <<"  " << m_max_nForce <<endl;
    cout << "m_shForce= " <<m_shForce <<"  " << m_max_shForce <<endl;
    cout << "m_tMoment= " <<m_tMoment <<"  " << m_max_tMoment <<endl;
    cout << "m_bMoment= " <<m_bMoment <<"  " << m_max_bMoment <<endl; }
    */



//    if (m_nForce<0.0) shear_tensile = 1 ;

//  cout <<  (m_p1->getOldPos() + m_p2->getOldPos())*0.5<<"  "<<   0.5*m_nForce*m_nForce/m_kr
//                      +0.5*m_shForce*m_shForce/m_ks
//                      +0.5*m_bMoment*m_bMoment/m_kb << "  " << shear_tensile <<  endl ;


  } else {res= false; }
  return res ;
}




Vec3 CRotThermBondedInteraction::getInitialCentrePtDiff() const
{
  return m_p2->getInitPos() - m_p1->getInitPos();
}

Vec3 CRotThermBondedInteraction::getCentrePtDiff() const
{
  return m_p2->getPos() - m_p1->getPos();
}

Vec3 CRotThermBondedInteraction::getInitialMidPoint() const
{
  const Vec3   initialDiff     = getInitialCentrePtDiff();
  const double initialDiffNorm = initialDiff.norm();
  return (m_p1->getRad() + (initialDiffNorm - m_p2->getRad()))/(2.0*initialDiffNorm)*initialDiff;
}

double CRotThermBondedInteraction::getCriterion() const
{
  return m_nForce/m_max_nForce   +
    m_shForce/m_max_shForce +
    m_tMoment/m_max_tMoment +
    m_bMoment/m_max_bMoment;
}

Vec3 CRotThermBondedInteraction::getShearDiff() const
{
  const Vec3 initialMidPt = getInitialMidPoint();
  const Vec3 p1Pt = m_p1->getPos() + (m_p1->getQuat().to_matrix()*(initialMidPt));
  const Vec3 p2Pt = m_p2->getPos() + (m_p2->getQuat().to_matrix()*(initialMidPt-m_p2->getInitPos()));
  const Vec3 ptDiff = p2Pt-p1Pt;
  const Vec3 rotatedInitialDiffDirection  = p1Pt - m_p1->getPos();  
  return (ptDiff - (((ptDiff*rotatedInitialDiffDirection)/rotatedInitialDiffDirection.norm2())*rotatedInitialDiffDirection));
}

#if 0
void CRotBondedInteraction::calcForces()
{
  // Calculate linear elastic force, no moment
  const Vec3 currentDiff = m_p2->getPos() - m_p1->getPos();
  const double currentDiffNorm = currentDiff.norm();
  const Vec3 initialDiff = getInitialCentrePtDiff();
  const Vec3 p2LinearForce = ((m_kr * (initialDiff.norm() - currentDiffNorm))/currentDiffNorm) * currentDiff;
  
  m_p2->applyForce( p2LinearForce, m_p2->getPos());
  m_p1->applyForce(-p2LinearForce, m_p1->getPos());
  m_nForce = p2LinearForce.norm();
  
  // Calculate the shearing forces, linear and moment contributions.
  // Need to take into account relative orientations of particles.
  const Vec3 shearDiff = getShearDiff();
}

#else

/*
const double HALF_SQRT_2 = 0.5*sqrt(2.0);
void CRotBondedInteraction::calcForces()
{
  double s_fai=0.0, c_fai=0.0, c_pasi2=0.0, s_pasi2=0.0;
  double sita=0.0,  pasi=0.0;

  const Matrix3 mat0 = (m_p1->getQuat()).to_matrix();
  const Vec3 rbp  = m_p2->getPos() - m_p1->getPos();
  const double rbpNorm = rbp.norm();
  const Vec3 rb = mat0*(m_p2->getPos() - m_p1->getPos());
  const double rbNorm = rb.norm();

  const Vec3 r_0 =  m_p2->getInitPos() - m_p1->getInitPos();
  const double r_0Norm = r_0.norm();
  const Vec3 delta_r = (rbNorm - r_0Norm)*rb/rbNorm;
  const Vec3 Fr = m_kr*delta_r ;

  const double gama_cos = dot(rb, r_0)/(rbNorm*r_0Norm) ; 
  double gama = 0.0;
  if (gama_cos < 1.0 && gama_cos > -1.0 )
  {
    gama = acos(gama_cos);
  }

  const Vec3 rb_0 = cross(rb, r_0);
  const Vec3 rb_b0 = cross(rb, rb_0);
  const double rb_b0Norm = rb_b0.norm();
  Vec3 Fst(0.0,0.0,0.0);
  if (rb_b0Norm != 0)
  {
    Fst = m_ks*r_0Norm*gama*rb_b0/rb_b0Norm;
  }

  Vec3 Mst(0.0,0.0,0.0);
  const double rb_0Norm = rb_0.norm();
  if (rb_0Norm != 0)
  {
    Mst = -Fst.norm()*rb_0/rb_0Norm;
  }
  
  const double m0 = HALF_SQRT_2*sqrt((rbNorm + rb.Z())/rbNorm);
  double m1       = -HALF_SQRT_2*sqrt((rbNorm - rb.Z())/rbNorm) ;
  double m2       = -m1;
  const double m3 = 0.0;
  const double rbX2PlusRbY2 = rb.X()*rb.X() + rb.Y()*rb.Y();
  if (rbX2PlusRbY2 != 0.0)
  {
    const double denomSqrt = sqrt(rbX2PlusRbY2);
    m1 = m1*rb.Y()/denomSqrt;
    m2 = m2*rb.X()/denomSqrt;
  }

  const Quaternion qm(m0, Vec3(m1,m2,m3));
  const Matrix3 mat   = qm.to_matrix();
  const Matrix3 matTrans(mat.trans());
  const Quaternion rp = (m_p1->getQuat()).inverse() * (m_p2->getQuat());
  const Quaternion r = qm.inverse() * rp * qm;
  const double temp0 = r.return_sca()*r.return_sca() +
                       r.return_vec().Z()*r.return_vec().Z();
  if(temp0 == 0.0)
  {
    pasi = 0.0;
  }
  else
  {
    const double sqrtTemp0 = sqrt(temp0);
    c_pasi2 = r.return_sca()/sqrtTemp0;
    s_pasi2 = r.return_vec().Z()/sqrtTemp0;
    pasi = 2.0 * calc_angle(s_pasi2,c_pasi2);
  }
  const Vec3 Mtp = Vec3(0.0,0.0, m_kt*pasi);
  const Vec3 Mtr = matTrans * Mtp;
  const double temp2 = r.return_vec().X()*r.return_vec().X() +
                       r.return_vec().Y()*r.return_vec().Y();
  const double temp3 = r.return_vec().X()*r.return_vec().Z() +
                       r.return_sca() * r.return_vec().Y();
  const double temp4 = r.return_vec().Y()*r.return_vec().Z() -
                       r.return_sca() * r.return_vec().X();

  Vec3 Mbr(0, 0, 0);
  Vec3 Fsr(0, 0, 0);
  Vec3 Msr(0, 0, 0);                       
  if (temp2 != 0.0)
  {
    const double temp1 = temp0 -  temp2;
    if (temp1 >= 1.0 || temp1 <= -1.0)
    {
      sita = 0.0 ;
    }
    else
    {
      sita = acos(temp1);
    }
    if (temp0 == 0.0) {
      const double sqrtDenom = sqrt(temp2);
      c_fai =  r.return_vec().Y()/sqrtDenom;
      s_fai = -r.return_vec().X()/sqrtDenom;
    } else {
      const double sqrtDenom = sqrt(temp0*temp2);
      c_fai = temp3/sqrtDenom;
      s_fai = temp4/sqrtDenom;
    }
    const Vec3 Mbp = Vec3( - 0.5*m_kb*sita*s_fai, 0.5*m_kb*sita*c_fai, 0.0);
    Mbr = matTrans * Mbp;
    const Vec3 Fsp = Vec3 ( - 0.5*m_ks*c_fai*sita*r_0Norm,
                            - 0.5*m_ks*s_fai*sita*r_0Norm,
                              0.0);
    Fsr = matTrans * Fsp;

    const Vec3 Msp = Vec3(   0.5*m_ks*s_fai*sita*r_0Norm,
                           - 0.5*m_ks*c_fai*sita*r_0Norm,
                             0.0  );
    Msr = matTrans * Msp;
  }

  const double eq_rad1 = m_p1->getRad()/(m_p1->getRad()+m_p2->getRad());
  const double eq_rad2 = m_p2->getRad()/(m_p1->getRad()+m_p2->getRad());

  const Matrix3 mat0Trans(mat0.trans());
  m_force  =  mat0Trans * (Fr + Fst + Fsr );
  m_moment =  mat0Trans * (Mbr + Mtr + Msr*eq_rad1*rbpNorm + Mst*eq_rad1*rbpNorm);

  const Vec3 moment2 = mat0Trans *( Msr*eq_rad2*rbpNorm + Mst*eq_rad2*rbpNorm - Mbr - Mtr);  

  m_shForce = (Fst + Fsr).norm();
  m_tMoment = Mtr.norm();
  m_bMoment = Mbr.norm();
  const double r0 = m_p1->getRad()+m_p2->getRad();
  const Vec3 D = m_p1->getPos() - m_p2->getPos();
  m_dist = sqrt(D*D);

  m_nForce = ((m_dist < r0) ? -1.0 : 1.0) * Fr.norm();

  Vec3 pos=m_p2->getPos()+(m_p2->getRad()/(m_p1->getRad()+m_p2->getRad()))*D;

  m_p2->applyForce(-1.0*m_force,pos);
  m_p1->applyForce(m_force,pos);
  m_cpos=pos;

  m_p1->applyMoment(m_moment) ;
  m_p2->applyMoment(moment2) ;
}
*/

void CRotThermBondedInteraction::calcForces()
{
   double  s_fai=0.0, c_fai=0.0, c_pasi2=0.0, s_pasi2=0.0;
   double sita=0.0,  pasi=0.0 ;
   //double temp_m1=0.0, temp_m2=0.0 ;
   double m0=0.0,m1=0.0,m2=0.0,m3=0.0;
   Vec3  Fst, Mst ;
   double temp0=0.0, temp1=0.0,temp2=0.0,temp3=0.0, temp4=0.0;
   Vec3 Mtr, Mbr, Fsr, Msr ;

   Mbr = Vec3 (0.0,0.0,0.0) ;
   Fsr = Vec3 (0.0,0.0,0.0) ;
   Msr = Vec3 (0.0,0.0,0.0) ;
  

//zjr 
//      cout << "wyc  in CRotBI   " << m_p1->getID() <<"--  " <<m_p2->getID() << endl;
//  cout << m_p1->getPos() << "  pos " << m_p2->getPos() << endl;
 // cout << m_p1->getQuat().return_sca() << "  "<< m_p1->getQuat().return_vec() <<" 1 quat " << endl;
 // cout << m_p2->getQuat().return_sca() << "  "<< m_p2->getQuat().return_vec() <<" 2 quat " << endl;

 //cout << " k=  " << m_kr << "  " << m_ks<<"  " << m_kt << "  " << m_kb << endl; 

  Matrix3 mat0 = (m_p1->getQuat()).to_matrix() ;

  
  //cout << " wyc  mat0 =  " << mat0 << endl;

   //Vec3 rbp  = m_p2->getPos() - m_p1->getPos() ;
   Vec3 rb = mat0*(m_p2->getPos() - m_p1->getPos() ) ;


   Vec3 r_0 =  m_p2->getInitPos() - m_p1->getInitPos() ;
   Vec3 delta_r = ( rb.norm()-r_0.norm() )*rb/rb.norm() ;
  // cout << " rb= "<< rb.norm()-r_0.norm() << endl;


   Vec3 Fr = m_kr*delta_r ;


   double gama_cos = dot(rb, r_0)/(rb.norm() *r_0.norm()) ;

   double gama;
    if (gama_cos >=1.0|| gama_cos<=-1.0 )   gama = 0.0 ;
    else   gama = acos(gama_cos);


   Vec3 rb_0 = cross(rb, r_0) ;
   Vec3 rb_b0 = cross(rb, rb_0) ;  //cout <<" 0 ?  "  << rb_b0.norm() << endl;
   if (rb_b0.norm() ==0){  Fst = Vec3(0.0,0.0,0.0) ; }
//   else {  Fst = m_ks*rb.norm()*gama* rb_b0/rb_b0.norm() ; }
// change oct 6/2005
     else {  Fst = m_ks*r_0.norm()*gama* rb_b0/rb_b0.norm() ; }     

  if (rb_0.norm() ==0) { Mst = Vec3(0.0,0.0,0.0) ;}
   else {

           Mst = -Fst.norm()*rb.norm()*rb_0/rb_0.norm()  ;
         }
   m0 = sqrt( (rb.norm()+ rb.Z() )/rb.norm() );
   if( rb.X()*rb.X() + rb.Y()*rb.Y() ==0.0 ) {
           m1 = 0.0 ; //-sqrt( (rb.norm()-rb.Z())/rb.norm() ) ;
           m2 = 0.0 ; // sqrt( (rb.norm()-rb.Z())/rb.norm() ) ;

   }else { //temp_m1 = rb.Y()/sqrt(rb.X()*rb.X()+rb.Y()*rb.Y());
           //temp_m2 = rb.X()/sqrt(rb.X()*rb.X()+rb.Y()*rb.Y());
          m1 = -sqrt( (rb.norm()-rb.Z())/rb.norm() )*rb.Y()/sqrt(rb.X()*rb.X()+rb.Y()*rb.Y());
          m2 =  sqrt( (rb.norm()-rb.Z())/rb.norm() )*rb.X()/sqrt(rb.X()*rb.X()+rb.Y()*rb.Y());
   }

   m3 =0.0 ;


   Quaternion qm(m0, Vec3(m1,m2,m3) ) ;
   Matrix3 mat   = (qm.to_matrix())*0.5 ;

 Quaternion rp = (m_p1->getQuat()).inverse() * (m_p2->getQuat())  ;

  double r_0_tmp = rp.return_sca() ;
  double r_1_tmp = rp.return_vec().X()*(m0*m0+m1*m1-m2*m2)
                  +2.0*m1*m2*rp.return_vec().Y()
                  -2.0*m0*m2*rp.return_vec().Z() ;
  double r_2_tmp = rp.return_vec().Y()*(m0*m0-m1*m1+m2*m2)
                  +2.0*m1*m2*rp.return_vec().X()
                  +2.0*m0*m1*rp.return_vec().Z() ;
  double r_3_tmp = rp.return_vec().Z()*(m0*m0-m1*m1-m2*m2)
                  +2.0*m0*m2*rp.return_vec().X()
                  -2.0*m0*m1*rp.return_vec().Y() ;
  Quaternion r(r_0_tmp, Vec3( r_1_tmp, r_2_tmp, r_3_tmp)*0.5) ;


   temp0 = r.return_sca()*r.return_sca() +
           r.return_vec().Z()*r.return_vec().Z();
   if( temp0 == 0.0 ) {  pasi = 0.0;}
   else {
         c_pasi2 = r.return_sca()/sqrt(temp0);
         s_pasi2 = r.return_vec().Z()/ sqrt(temp0);
         pasi = 2.0 * calc_angle(s_pasi2, c_pasi2);
   }
   Vec3 Mtp = Vec3(0.0,0.0, m_kt*pasi) ;
   Mtr = mat.trans() * Mtp ;
   temp2 = r.return_vec().X()*r.return_vec().X() +
           r.return_vec().Y()*r.return_vec().Y()    ;
   temp3 = r.return_vec().X()*r.return_vec().Z() +
           r.return_sca() * r.return_vec().Y()      ;
   temp4 = r.return_vec().Y()*r.return_vec().Z() -
           r.return_sca() * r.return_vec().X()      ;


   if ( temp2 != 0.0) {
      temp1 =  r.return_sca()*r.return_sca()
             + r.return_vec().Z()*r.return_vec().Z()
             - r.return_vec().X()*r.return_vec().X()
             - r.return_vec().Y()*r.return_vec().Y() ;
      if (temp1 >=1.0|| temp1<=-1.0 )  { sita = 0.0 ;}
      else { 
 
               sita = acos(temp1);}
      if( temp0  == 0.0) {
         c_fai = r.return_vec().Y()/sqrt(temp2);
         s_fai = -r.return_vec().X()/sqrt(temp2);


      }else {
         c_fai = temp3/sqrt(temp0 * temp2) ;
         s_fai = temp4/sqrt(temp0 * temp2) ;
      }
      Vec3 Mbp = Vec3( - 0.5*m_kb*sita*s_fai, 0.5*m_kb*sita*c_fai, 0.0) ;
      Mbr = mat.trans() * Mbp  ;
      Vec3 Fsp = Vec3 ( - m_ks*c_fai*sita*m_p2->getRad() ,
                        - m_ks*s_fai*sita*m_p2->getRad() ,
                          0.0) ;
      Fsr = mat.trans() * Fsp  ;


       Vec3 Msp = Vec3(   m_ks*s_fai*sita*m_p2->getRad()*rb.norm() ,
                        - m_ks*c_fai*sita*m_p2->getRad()*rb.norm()  ,
                          0.0  );


      Msr = mat.trans() * Msp   ;
   }
   double eq_rad1 = m_p1->getRad()/(m_p1->getRad()+m_p2->getRad()) ;
   double eq_rad2 = m_p2->getRad()/(m_p1->getRad()+m_p2->getRad()) ;


    m_force  =  mat0.trans() *(Fr + Fst + Fsr ) ;

//   cout << "fr= " << Fr <<" Fst= " << Fst <<" Fsr= " <<Fsr<<endl;
 
//zjr   
//     cout << "wyc m_force= " << m_force << endl;

//    m_moment =  mat0.trans() *(Mbr + Mtr + Msr*eq_rad1*rbp.norm() + Mst*eq_rad1*rbp.norm()) ;
// oct 6 /2005 changed

    m_moment =  mat0.trans() *(Mbr + Mtr + Msr*eq_rad1 + Mst*eq_rad1) ;

//zjr  
//      cout << "wyc  m_moment=" <<m_moment << endl;
//zjr  
//     cout << "bending angle=  " << sita<< "  twisting angle=  " << pasi << endl;

//2d case , temp


//  if( m_force.Z()!=0.0) m_force = Vec3( m_force.X(), m_force.Y(), 0.0) ;

//  if(m_moment.X()!=0.0) m_moment = Vec3(0.0, m_moment.Y(), m_moment.Z()) ;
//  if(m_moment.Y()!=0.0) m_moment = Vec3(m_moment.X(),0.0, m_moment.Z()) ;

 

   Vec3 moment2 = mat0.trans() *( Msr*eq_rad2 + Mst*eq_rad2 - Mbr - Mtr) ;
 
// 2d case, temp

//  if(moment2.X()!=0.0) moment2 = Vec3(0.0, moment2.Y(), moment2.Z()) ;
//  if(moment2.Y()!=0.0) moment2 = Vec3(moment2.X(),0.0, moment2.Z()) ;



    m_shForce = (Fst + Fsr).norm();
   m_tMoment = Mtr.norm();
   m_bMoment = Mbr.norm();
   double Eq_r0= m_p1->getRad()+m_p2->getRad();
   Vec3 D=m_p1->getPos()-m_p2->getPos();
   m_dist=sqrt(D*D);
   if(m_dist<Eq_r0) {m_nForce = -1.0*Fr.norm();  }
   else {m_nForce = Fr.norm();}
   
   D=m_p1->getPos()-m_p2->getPos();

   Vec3 pos=m_p2->getPos()+(m_p2->getRad()/(m_p1->getRad()+m_p2->getRad()))*D;

   m_p2->applyForce(-1.0*m_force,pos);
   m_p1->applyForce(m_force,pos);
   m_cpos=pos;

   m_p1->applyMoment(m_moment) ;
   m_p2->applyMoment(moment2) ;

}
#endif

Vec3 CRotThermBondedInteraction::getForce() const
{
  return m_force;
}

double CRotThermBondedInteraction::getPotentialEnergy() const
{
  double pe_r,pe_s,pe_t,pe_b;

  pe_r=(m_kr!=0.0) ? 0.5*(m_nForce*m_nForce)/m_kr : 0.0;
  pe_s=(m_ks!=0.0) ? 0.5*(m_shForce*m_shForce)/m_ks : 0.0;
  pe_t=(m_kt!=0.0) ? 0.5*(m_tMoment*m_tMoment)/m_kt : 0.0;
  pe_b=(m_kb!=0.0) ? 0.5*(m_bMoment*m_bMoment)/m_kb : 0.0;
  return pe_r+pe_s+pe_t+pe_b;
}

double CRotThermBondedInteraction::getNormalPotentialEnergy() const
{
  return (m_kr!=0.0) ? 0.5*(m_nForce*m_nForce)/m_kr : 0.0;
}

double CRotThermBondedInteraction::getShearPotentialEnergy() const
{
  return (m_ks!=0.0) ? 0.5*(m_shForce*m_shForce)/m_ks : 0.0;
}
 
double CRotThermBondedInteraction::getTwistPotentialEnergy() const
{
  return (m_kt!=0.0) ? 0.5*(m_tMoment*m_tMoment)/m_kt : 0.0;
}

double CRotThermBondedInteraction::getBendPotentialEnergy() const
{
  return (m_kb!=0.0) ? 0.5*(m_bMoment*m_bMoment)/m_kb : 0.0;
}

Vec3 CRotThermBondedInteraction::getBondedVector1() const
{

  Vec3 D =  m_p2->getPos()-m_p1->getPos() ;
  Vec3 pos = m_p1->getPos() + D*m_p1->getRad()/(m_p1->getRad()+m_p2->getRad());

  return m_p1->getPos()-pos ;

}

Vec3 CRotThermBondedInteraction::getBondedVector2() const
{

  Vec3 D =  m_p2->getPos()-m_p1->getPos() ;
  Vec3 pos = m_p1->getPos() + D*m_p1->getRad()/(m_p1->getRad()+m_p2->getRad());

  return m_p2->getPos()-pos ; 
}

/*!
  Get the particle member function which returns a scalar field of a given name.

  \param name the name of the field 
*/
CRotThermBondedInteraction::ScalarFieldFunction CRotThermBondedInteraction::getScalarFieldFunction(const string& name)
{
  CRotThermBondedInteraction::ScalarFieldFunction sf;
                                                                                
  if (name=="potential_energy"){
    sf=&CRotThermBondedInteraction::getPotentialEnergy;
  } else if (name=="e_pot_normal"){
    sf=&CRotThermBondedInteraction::getNormalPotentialEnergy;
  } else if (name=="e_pot_shear"){
    sf=&CRotThermBondedInteraction::getShearPotentialEnergy;
  } else if (name=="e_pot_twist"){
    sf=&CRotThermBondedInteraction::getTwistPotentialEnergy;
  } else if (name=="e_pot_bend"){
    sf=&CRotThermBondedInteraction::getBendPotentialEnergy;
  } else if (name=="count"){
    sf=&CRotThermBondedInteraction::Count;
  } else if (name=="breaking_criterion"){
    sf=&CRotThermBondedInteraction::getCriterion;
  } else {
    sf=NULL;
    cerr << "ERROR - invalid name for interaction scalar access function" << endl;
  }
                                                                                
  return sf;
}

/*!
  Get the particle member function which returns a vector field of a given name.

  \param name the name of the field 
*/
CRotThermBondedInteraction::VectorFieldFunction CRotThermBondedInteraction::getVectorFieldFunction(const string& name)
{
  CRotThermBondedInteraction::VectorFieldFunction vf;
                                                                                
  if (name=="force"){
    vf=&CRotThermBondedInteraction::getForce;
  } else {
    vf=NULL;
    cerr << "ERROR - invalid name for interaction vector access function" << endl;
  }
                                                                                
  return vf;
}

/*!
  Get the particle member function which returns a checked scalar field of a given name.

  \param name the name of the field
*/
CRotThermBondedInteraction::CheckedScalarFieldFunction CRotThermBondedInteraction::getCheckedScalarFieldFunction(const string& name)
{
  CRotThermBondedInteraction::CheckedScalarFieldFunction sf;

  sf=NULL;
  cerr << "ERROR - invalid name for interaction scalar access function" << endl;

  return sf;
}

void CRotThermBondedInteraction::saveCheckPointData(std::ostream &oStream)
{
  throw std::runtime_error(
    "Checkpointing not implemented for CRotThermBondedInteraction."
  );
}

void CRotThermBondedInteraction::loadCheckPointData(std::istream &iStream)
{
  throw std::runtime_error("CRotThermBondedInteraction::loadCheckPointData not implemented.");
}

void CRotThermBondedInteraction::calcHeatTrans()
{

  double  eta = 1.5 ;

  Vec3 D=m_p1->getPos()-m_p2->getPos();
  double Rij2 = D.norm2();
  double d_temp = m_p2->getTemperature() -  m_p1->getTemperature() ; 
  double heatij = eta*m_diffusivity*d_temp/Rij2 ;
  m_p1->applyHeatTrans(heatij) ;
  m_p2->applyHeatTrans(-heatij) ;
}

/*!
  Pack this object into a TML packed message

  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::pack<CRotThermBondedInteraction>(const CRotThermBondedInteraction& I)
{
  append(I.m_id[0]);
  append(I.m_id[1]);
  append(static_cast<int>(I.m_init));
  append(I.m_kr);
  append(I.m_ks);
  append(I.m_kb);
  append(I.m_kt);
  append(I.m_max_nForce);
  append(I.m_max_shForce);
  append(I.m_max_tMoment);
  append(I.m_max_bMoment);
  append(I.m_diffusivity);
}

/*!
  Unpack a CBondedInteraction from a TML packed message
                                                                                
  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::unpack<CRotThermBondedInteraction>(CRotThermBondedInteraction& I)
{
  I.m_id.clear();
  I.m_id.push_back(pop_int());
  I.m_id.push_back(pop_int());
  I.m_init = static_cast<bool>(pop_int());
  I.m_kr = pop_double();
  I.m_ks = pop_double();
  I.m_kb = pop_double();
  I.m_kt = pop_double();
  I.m_max_nForce = pop_double();
  I.m_max_shForce = pop_double();
  I.m_max_tMoment = pop_double();
  I.m_max_bMoment = pop_double();
  I.m_diffusivity = pop_double();
}

/*!
  save restart data to ostream

  \param oStream the output stream
*/
void CRotThermBondedInteraction::saveRestartData(std::ostream &oStream)
{
  oStream << m_id[0] << " ";
  oStream << m_id[1] << " ";
  oStream << m_init  << " ";
  oStream << m_kr << " ";
  oStream << m_ks << " ";
  oStream << m_kb << " ";
  oStream << m_kt << " ";
  oStream << m_max_nForce << " ";
  oStream << m_max_shForce << " ";
  oStream << m_max_tMoment << " ";
  oStream << m_max_bMoment << " ";
  oStream << m_diffusivity << " ";
  oStream << getTag() ;
}

/*!
  load restart data from stream

  \param iStream the input stream
*/
void CRotThermBondedInteraction::loadRestartData(std::istream &iStream)
{
  int tag;
  iStream >> m_id[0];
  iStream >> m_id[1];
  iStream >> m_init ;
  iStream >> m_kr;
  iStream >> m_ks;
  iStream >> m_kb;
  iStream >> m_kt;
  iStream >> m_max_nForce;
  iStream >> m_max_shForce;
  iStream >> m_max_tMoment;
  iStream >> m_max_bMoment;
  iStream >> m_diffusivity;
  iStream >> tag;
  setTag(tag);
}

ostream& operator<<(ostream& ost,const CRotThermBondedInteraction& BI)
{
  ost << "[" << BI.m_p1->getID() << " - " << BI.m_p2->getID() << "]";
  return ost;
}
