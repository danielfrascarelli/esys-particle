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

#include "Model/ARotBondedInteraction.h"
#include "Model/BondedInteractionCpData.h"
#include "Foundation/console.h"

// --- TML includes ---
#include "tml/message/packed_message_interface.h"

#include <stdexcept>

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


/*!
    Default constructor for basic rotational bonded interaction parameters.
    Sets all stiffnesses to 0.0
*/
ARotBondedIGP::ARotBondedIGP()
 : AIGParam(),
   kr(0.0),
   ks(0.0),
   kt(0.0),
   kb(0.0),
   tag(0),
   scaling(true)
{}

/*!
    Constructor for basic rotational bonded interaction parameters from
    individual stiffness parameters (Wang et al 2006)
*/
ARotBondedIGP::ARotBondedIGP(
  const  std::string &name,
  double kr,
  double ks,
  double kt,
  double kb,
  int    tag,
  bool   scaling
)
 : AIGParam(name),
   kr(kr),
   ks(ks),
   kt(kt),
   kb(kb),
   tag(tag),
   scaling(scaling)
{}

/*!
    Constructor for basic rotational bonded interaction parameters from
    brittle beam parameters (stiffness only, not strength).
*/
ARotBondedIGP::ARotBondedIGP(
  const  std::string &name,
  double youngsModulus,
  double poissonsRatio,
  int    tag
)
 : AIGParam(name),
   tag(tag),
   scaling(true)
{
   double shearModulus = youngsModulus / (2.0*(1.+poissonsRatio));

   kr = M_PI*youngsModulus;
   ks = M_PI*shearModulus;
   kb = M_PI*youngsModulus;
   kt = M_PI*shearModulus;
}

// --------- END PARAMETER CLASS ----------


/*!
    Default constructor for base of rotational bonded interactions.
    Just sets all parameters to 0.0 
*/
ARotBondedInteraction::ARotBondedInteraction():ARotPairInteraction()
{
  m_kr = 0.0 ;
  m_ks = 0.0 ;
  m_kb = 0.0 ;
  m_kt = 0.0 ;

  m_nForce  = 0.0 ;
  m_shForce = 0.0;
  m_tMoment = 0.0;
  m_bMoment = 0.0;
  m_moment  = Vec3(0.0,0.0,0.0);
  m_force   = Vec3(0.0,0.0,0.0);
  
  m_tag = 0;
  m_scaling = true;
}

/*!
    Constructor for base of rotational bonded interactions using two particle pointers
    and a parameter set contained in a ARotBondedIGP.
*/
ARotBondedInteraction::ARotBondedInteraction(CRotParticle* p1,CRotParticle* p2,const ARotBondedIGP& param):ARotPairInteraction(p1,p2)
{
  m_nForce  = 0.0;
  m_shForce = 0.0;
  m_tMoment = 0.0;
  m_bMoment = 0.0;

  m_D = p2->getPos()-p1->getPos();
  m_dist=sqrt(m_D*m_D);
  m_r0=m_dist;

  double effR=1.0;
  double effL=1.0;
  double effA=1.0;
  double momentJ=1.0;
  double momentI=1.0;
  m_scaling = param.scaling;
 
  // scale elastic param
  // assume Rmean-scaling -> add option to do otherwise to derived class
  if (m_scaling == true) {
    if(!CParticle::getDo2dCalculations()){
      effR=0.5*m_dist;
      effL = m_dist;
      effA = effR * effR;
      momentJ = 0.5 * effR * effR * effR * effR;
      momentI = 0.5 * momentJ;
    }
  }

  m_kr =  param.kr * effA / effL;
  m_ks =  param.ks * effA / effL;
  m_kb =  param.kb * momentI / effL;
  m_kt =  param.kt * momentJ / effL;
  
  m_force  = Vec3(0.0,0.0,0.0);
  m_moment = Vec3(0.0,0.0,0.0) ;

  m_tag = param.tag;
}

/*!
    Default destructor - does nothing
*/
ARotBondedInteraction::~ARotBondedInteraction()
{}

/*!
    get bond tag
*/
int ARotBondedInteraction::getTag() const
{
  return m_tag;
}

/*!
    set bond tag
*/
void ARotBondedInteraction::setTag(int tag)
{
  m_tag = tag;
}

/*!
    Get _initial_ vector between particle centers (i.e. at contruction time)
*/
Vec3 ARotBondedInteraction::getInitialCentrePtDiff() const
{
  return m_p2->getInitPos() - m_p1->getInitPos();
}

/*!
    Get _current_ vector between particle centers
*/
Vec3 ARotBondedInteraction::getCentrePtDiff() const
{
  return m_p2->getPos() - m_p1->getPos();
}

/*!
    Get location of the midpoint between the particles at construction time
*/
Vec3 ARotBondedInteraction::getInitialMidPoint() const
{
  const Vec3   initialDiff  = getInitialCentrePtDiff();
  const Vec3 normalisedDiff = initialDiff*(1.0/initialDiff.norm());
  return 0.5*((m_p1->getInitPos() + normalisedDiff*m_p1->getRad()) + (m_p2->getInitPos() - normalisedDiff*m_p2->getRad()));
}

/*!
    Get shear force on particle 2
*/
Vec3 ARotBondedInteraction::getP2ShearForcePt() const
{
  return m_p2->getPos() + (m_p2->getQuat().to_matrix()*(getInitialMidPoint()-m_p2->getInitPos()));
}

/*!
    Get shear force on particle 1
*/
Vec3 ARotBondedInteraction::getP1ShearForcePt() const
{
  return m_p1->getPos() + (m_p1->getQuat().to_matrix()*(getInitialMidPoint()-m_p1->getInitPos()));
}

/*!
    ?
*/
Vec3 ARotBondedInteraction::getShearDiff() const
{
  const Vec3 p1Pt = getP1ShearForcePt();
  const Vec3 p2Pt = getP2ShearForcePt();
  const Vec3 ptDiff = p2Pt-p1Pt;
  const Vec3 rotatedInitialDiffDirection  = p1Pt - m_p1->getPos(); //getCentrePtDiff();
  return (ptDiff - (((ptDiff*rotatedInitialDiffDirection)/rotatedInitialDiffDirection.norm2())*rotatedInitialDiffDirection));
}

/*!
    get _current_ midpoint between particles
*/
Vec3 ARotBondedInteraction::getContactPoint() const
{
  const Vec3 centrePtDiff = getCentrePtDiff();
  const Vec3 normalisedDiff = centrePtDiff*(1.0/centrePtDiff.norm());
  return 0.5*((m_p1->getPos() + normalisedDiff*m_p1->getRad()) + (m_p2->getPos() - normalisedDiff*m_p2->getRad()));
}

#if 0

Vec3 getTangent(const Vec3 &r, const Vec3 &f)
{
  return (f.norm2() > 0) ? (f - ((r*f)/f.norm2())*r) : Vec3::ZERO;
}

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
  const Vec3 p2ShearForcePt = getP2ShearForcePt();
  const Vec3 p1ShearForcePt = getP1ShearForcePt();
  const Vec3 p2ShearForce = -m_ks * (p2ShearForcePt - p1ShearForcePt); //(m_ks * shearDiff);
  m_p2->applyMoment(cross(p2ShearForcePt - m_p2->getPos(), getTangent(p2ShearForcePt - m_p2->getPos(), p2ShearForce)));
  m_p1->applyMoment(cross(p1ShearForcePt - m_p1->getPos(), getTangent(p1ShearForcePt - m_p1->getPos(), -p2ShearForce)));
  m_p2->applyForce( p2ShearForce, m_p2->getPos());
  m_p1->applyForce(-p2ShearForce, m_p1->getPos());
  m_shForce = p2ShearForce.norm();
}

#else

const double HALF_SQRT_2 = 0.5*sqrt(2.0); // use constexpr?

/*!
    Calculate forces
*/
void ARotBondedInteraction::calcForces()
{
  double s_fai=0.0, c_fai=0.0, c_pasi2=0.0, s_pasi2=0.0;
  double sita=0.0,  pasi=0.0;

  const Matrix3 mat0 = (m_p1->getQuat()).to_matrix();
  const Vec3 rbp  = m_p2->getPos() - m_p1->getPos();  // rbp <-> r_f 
  const double rbpNorm = rbp.norm();
  const Vec3 rb = mat0*(m_p2->getPos() - m_p1->getPos()); // eq. 11, rb <-> r_c
  const double rbNorm = rb.norm();

  const double r_0Norm = m_p1->getRad()+m_p2->getRad();  
  //const double r_0Norm = m_r0;
  const Vec3 delta_r = (rbNorm - r_0Norm)*rb/rbNorm; // eq. 9, part of eq. 12
  const Vec3 Fr = m_kr*delta_r ; // rest of eq. 12

  const double gama_cos = dot(rb, m_D)/(rbNorm*m_D.norm()) ; 
  double gama = 0.0;
  if (gama_cos < 1.0 && gama_cos > -1.0 )
  {
    gama = acos(gama_cos);
  }

  const Vec3 rb_0 = cross(rb, m_D);
  const Vec3 rb_b0 = cross(rb, rb_0);
  const double rb_b0Norm = rb_b0.norm();
  Vec3 Fst(0.0,0.0,0.0);
  if (rb_b0Norm != 0)
  {
    Fst = m_ks*r_0Norm*gama*rb_b0/rb_b0Norm; // eq. 14
  }

  Vec3 Mst(0.0,0.0,0.0);
  const double rb_0Norm = rb_0.norm();
  if (rb_0Norm != 0)
  {
    Mst = -Fst.norm()*rb_0/rb_0Norm; // eq. 15
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
    const Vec3 Mbp = Vec3( - 0.5*m_kb*sita*s_fai, 0.5*m_kb*sita*c_fai, 0.0); // eq 16.1
    Mbr = matTrans * Mbp; // eq. 17.1
    const Vec3 Fsp = Vec3 ( - 0.5*m_ks*c_fai*sita*r_0Norm,
                            - 0.5*m_ks*s_fai*sita*r_0Norm,
                0.0); // eq. 16.3
    Fsr = matTrans * Fsp; // eq. 17.3

    const Vec3 Msp = Vec3(   0.5*m_ks*s_fai*sita*r_0Norm,
                           - 0.5*m_ks*c_fai*sita*r_0Norm,
                             0.0  ); 
    // above is diff to paper: r instead r^2 here, other r moved to eq 19 below
    Msr = matTrans * Msp;
  }

  const double eq_rad1 = m_p1->getRad()/(m_p1->getRad()+m_p2->getRad());
  const double eq_rad2 = m_p2->getRad()/(m_p1->getRad()+m_p2->getRad());

  const Matrix3 mat0Trans(mat0.trans());
  m_force  =  mat0Trans * (Fr + Fst + Fsr );
  m_moment =  mat0Trans * (Mbr + Mtr + Msr*eq_rad1*rbpNorm + Mst*eq_rad1*rbpNorm); // eq. 19 for particle 1 

  const Vec3 moment2 = mat0Trans *( Msr*eq_rad2*rbpNorm + Mst*eq_rad2*rbpNorm - Mbr - Mtr);  // eq. 19 for particle 2

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

#endif

// --- FIELD FUNCTIONS ---
/*!
    get stored elastic energy
*/
double ARotBondedInteraction::getPotentialEnergy() const
{
  double pe_r,pe_s,pe_t,pe_b;

  pe_r=(m_kr!=0.0) ? 0.5*(m_nForce*m_nForce)/m_kr : 0.0;
  pe_s=(m_ks!=0.0) ? 0.5*(m_shForce*m_shForce)/m_ks : 0.0;
  pe_t=(m_kt!=0.0) ? 0.5*(m_tMoment*m_tMoment)/m_kt : 0.0;
  pe_b=(m_kb!=0.0) ? 0.5*(m_bMoment*m_bMoment)/m_kb : 0.0;
  return pe_r+pe_s+pe_t+pe_b;
}

double ARotBondedInteraction::getNormalPotentialEnergy() const
{
  return (m_kr!=0.0) ? 0.5*(m_nForce*m_nForce)/m_kr : 0.0;
}

double ARotBondedInteraction::getShearPotentialEnergy() const
{
  return (m_ks!=0.0) ? 0.5*(m_shForce*m_shForce)/m_ks : 0.0;
}
 
double ARotBondedInteraction::getTwistPotentialEnergy() const
{
  return (m_kt!=0.0) ? 0.5*(m_tMoment*m_tMoment)/m_kt : 0.0;
}

double ARotBondedInteraction::getBendPotentialEnergy() const
{
  return (m_kb!=0.0) ? 0.5*(m_bMoment*m_bMoment)/m_kb : 0.0;
}

Vec3 ARotBondedInteraction::getForce() const
{
  return m_force;
}

Vec3 ARotBondedInteraction::getNormalForce() const
{
  Vec3 dr = (m_p1->getPos() - m_p2->getPos()).unit();
  Vec3 normalForce = (m_force*dr)*dr;
  return normalForce;
}

Vec3 ARotBondedInteraction::getTangentialForce() const
{
  return (m_force - getNormalForce());
}

/*!
  Get the particle member function which returns a scalar field of a given name.

  \param name the name of the field 
*/
ARotBondedInteraction::ScalarFieldFunction ARotBondedInteraction::getScalarFieldFunction(const string& name)
{
  ARotBondedInteraction::ScalarFieldFunction sf;
                                                                                
  if (name=="potential_energy"){
    sf=&ARotBondedInteraction::getPotentialEnergy;
  } else if (name=="e_pot_normal"){
    sf=&ARotBondedInteraction::getNormalPotentialEnergy;
  } else if (name=="e_pot_shear"){
    sf=&ARotBondedInteraction::getShearPotentialEnergy;
  } else if (name=="e_pot_twist"){
    sf=&ARotBondedInteraction::getTwistPotentialEnergy;
  } else if (name=="e_pot_bend"){
    sf=&ARotBondedInteraction::getBendPotentialEnergy;
  } else if (name=="count"){
    sf=&ARotBondedInteraction::Count;
  } else {
    sf=NULL;
    cerr << "ERROR - invalid name for interaction scalar  access function" << endl;
  }
                                                                                
  return sf;
}

/*!
  Get the particle member function which returns a vector field of a given name.

  \param name the name of the field 
*/
ARotBondedInteraction::VectorFieldFunction ARotBondedInteraction::getVectorFieldFunction(const string& name)
{
  ARotBondedInteraction::VectorFieldFunction vf;
                                                                                
  if (name=="force"){
    vf=&ARotBondedInteraction::getForce;
  } else if (name=="normal_force"){
    vf=&ARotBondedInteraction::getNormalForce;
  } else if (name=="tangential_force"){
    vf=&ARotBondedInteraction::getTangentialForce;
  } else {
    vf=NULL;
    cerr << "ERROR - invalid name for interaction vector  access function" << endl; 
  }
                                                                                
  return vf;
}

/*!
  Get the particle member function which returns a checked scalar field of a given name.

  \param name the name of the field
*/
ARotBondedInteraction::CheckedScalarFieldFunction ARotBondedInteraction::getCheckedScalarFieldFunction(const string& name)
{
  ARotBondedInteraction::CheckedScalarFieldFunction sf;

  sf=NULL;
  cerr << "ERROR - invalid name for interaction scalar  access function" << endl;

  return sf;
}

