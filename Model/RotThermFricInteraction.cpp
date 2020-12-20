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
#include "Model/RotThermFricInteraction.h"
#include "Foundation/console.h"

// --- TML includes ---
#include "tml/message/packed_message_interface.h"

extern double calc_angle(double, double);

CRotThermFrictionIGP::CRotThermFrictionIGP()
: k(0.0),
  mu_d(0.0),
  mu_s(0.0),
  k_s(0.0),
  dt(0.0),
  diffusivity(0.0)
//  ds(Vec3(0.0,0.0,0.0))
{
}

CRotThermFrictionIGP::CRotThermFrictionIGP(
  const std::string &name,
  double normalK,
  double muDynamic,
  double muStatic,
  double shearK,
  double thermalDiffusivity,
  double deltaT
) :
  AIGParam(name), 
  k(normalK),
  mu_d(muDynamic),
  mu_s(muStatic),
  k_s(shearK),
  dt(deltaT),
  diffusivity(thermalDiffusivity)
{
}

CRotThermFrictionInteraction::CRotThermFrictionInteraction():ARotThermPairInteraction() 
{
  m_k=0.0;
  m_mu_d=0.0;
  m_mu_s=0.0;
  m_r0=0.0;
  m_ks=0.0;
  m_dt=0.0;
  m_is_slipping=false;
  m_is_touching=false;
  m_Ffric = Vec3(0.0,0.0,0.0);
  m_E_diss = 0.0; 
  m_diffusivity = 0.0;
//  m_ds = Vec3(0.0,0.0,0.0); 

  m_p1=NULL;
  m_p2=NULL;
  m_id.push_back(-1);
  m_id.push_back(-1);
}

CRotThermFrictionInteraction::CRotThermFrictionInteraction(
  CRotThermParticle* p1,
  CRotThermParticle* p2,
  const CRotThermFrictionIGP& param
) :
  ARotThermPairInteraction(p1,p2)
{
  m_mu_d = param.mu_d;
  m_mu_s = param.mu_s;
  m_r0=p1->getRad()+p2->getRad();
  m_dt=param.dt;
  m_cpos=p1->getPos()+((p2->getPos()-p1->getPos())*p1->getRad()/m_r0);
  m_is_slipping=false;
  m_is_touching=false;
  m_Ffric = Vec3(0.0,0.0,0.0);
  m_E_diss = 0.0; 
  m_ds = Vec3(0.0,0.0,0.0);

//wyc added   22/02/2005
  double min_r;
  double ran_ratio;
  double ran_ratioH;
  if (m_p1->getRad() <= m_p2->getRad()) min_r = m_p1->getRad();
  else                                  min_r = m_p2->getRad();
//  double ran_ratio = 2.0*min_r/(m_p1->getRad()+m_p2->getRad());
  if(m_p1->getDo2dCalculations()) { // 2D
    ran_ratio = 2.0*min_r/(m_p1->getRad()+m_p2->getRad());
    ran_ratioH = 2.0*min_r*(m_p1->getRad()+m_p2->getRad());
  } else { //  3D
    ran_ratio = 2.0*min_r*min_r/(m_p1->getRad()+m_p2->getRad());
    ran_ratioH = 2.0*min_r*min_r*(m_p1->getRad()+m_p2->getRad());
  }

  m_k = ran_ratio*param.k;
  m_ks = ran_ratio*param.k_s;
  m_diffusivity = ran_ratioH*param.diffusivity;
}

CRotThermFrictionInteraction::~CRotThermFrictionInteraction()
{
}

/*!
  Calculate elastic and frictional forces.
*/
void CRotThermFrictionInteraction::calcForces()
{  //cout << "wyc in RotFric:: calcF " <<endl;
  Vec3 pos;
  Vec3 force;
  Vec3 dv, ds;
  Vec3  d_Ffric;
  // calculate distance
  Vec3 D=m_p2->getPos()-m_p1->getPos();
  double dist=D*D;
  double eq_dist=m_p1->getRad()+m_p2->getRad();
  // check if there is contact
  if(dist<(eq_dist*eq_dist)){ // contact -> calculate forces
    //--- elastic force ---
    dist=sqrt(dist);
    force=D*(m_k*(dist-eq_dist)/dist);
    m_normal_force=force;
    pos=m_p2->getPos()-(m_p2->getRad()/eq_dist)*D;

    //cout << " wyc p1, p2 =  " << m_p1->getID() << "---  " << m_p2->getID()<<endl;   // cout << " wyc normal  =  " <<  force  << endl;

    // apply elastic force
    m_p2->applyForce(-1.0*force,pos);
    m_p1->applyForce(force,pos);

    //--- frictional force ---
    const Vec3 vp1_trs = m_p1->getVel();
    const Vec3 vp2_trs = m_p2->getVel();
    const Vec3 vp1_rot = cross(m_p1->getAngVel_t(),pos-m_p1->getPos());
    const Vec3 vp2_rot = cross(m_p2->getAngVel_t(),pos-m_p2->getPos());
    Vec3 dv_trs_s;

    const Vec3 dv_trs = vp2_trs-vp1_trs;
    // tangential part
    dv_trs_s = dv_trs - (dot(dv_trs,D)/D.norm2())*D;
    const Vec3 ds = (vp2_rot-vp1_rot + dv_trs_s )*m_dt;
    m_ds = ds;
//    cout<<" normal force = " << force.norm() << endl;
//    cout<< "ds= " <<  ds <<" norm= "<<ds.norm() <<  endl; //" trs= " <<dv_trs_s <<"rot=  " << vp2_rot-vp1_rot <<endl;
//    due to  motion of 2 particles as a rigid body !
//    Matrix3 mat0 = (m_p1->getQuat()).to_matrix() ;
//    cout  << "mat0= " <<mat0 <<"m_p1->getAngVel_t()=  " <<m_p1->getAngVel_t() <<   endl;
    Vec3 rbp  = m_p2->getPos() - m_p1->getPos() ;
//    Vec3 rb = mat0*rbp ;
    Vec3 vbp  = m_p2->getVel() - m_p1->getVel() ;
//    Vec3 vb = mat0*vbp ;
    double rbp0 = rbp.norm() ;
//    Vec3 omiga_s = mat0.trans()*m_p1->getAngVel_t();
//28/02/2005            getAngVel_t() : already in space-fixed system

    Vec3 omiga_s = 0.5*( m_p1->getAngVel_t()+m_p2->getAngVel_t());
    Vec3 omiga_spin = dot(omiga_s,rbp)*rbp/(rbp0*rbp0);
    Vec3 omiga_m = cross(rbp,vbp)/(rbp0*rbp0);
//    cout <<" omiga_spin= " <<omiga_spin << " omiga_m= " <<omiga_m << endl;
    d_Ffric = m_dt*cross(omiga_spin + omiga_m, m_Ffric);
//    cout << "D= " << D << "  d_Ffric=  " << d_Ffric<< "  m_Ffric= " <<m_Ffric  <<endl;
    m_Ffric += d_Ffric;

    if (m_is_slipping==false) {
//      if (!m_is_touching) {m_Ffric = Vec3(0.0,0.0,0.0);} //cout <<"first touch"<<endl;}//first touch
      if ((m_Ffric+m_ks*ds).norm()>force.norm()*m_mu_s) { // tangential force greater than static friction -> dynamic
        // m_Ffric= m_mu_d*force.norm() *ds/ds.norm();
        m_Ffric = m_mu_d*force.norm()*(m_Ffric+m_ks*ds)/(m_Ffric+m_ks*ds).norm();
        m_force_deficit = Vec3(0.0,0.0,0.0);
        m_is_slipping=true;
        m_E_diss = fabs(m_Ffric * ds); // energy dissipated
//        cout <<" stick->slip:  "<< endl;
//        cout << m_Ffric<<"  norm()= " <<  m_Ffric.norm()<<endl;
      } else { // static friction or no frictional force
        m_Ffric += m_ks*ds;
        m_E_diss = 0.0; // no energy dissipated
//        cout <<" continue stick :  " << " m_Ffric "<< m_Ffric<<"  norm= "  << m_Ffric.norm() <<  endl;
      }
      //  cout <<  "   m_Ffric=" <<  m_Ffric  << endl;
    } else if (m_is_slipping==true) {
//02/Nov  2005         it is hard to judge whether ds=0.0
//  criterion (1)
//      if(ds.norm()> 1.0e-8  ) {
//  criterion (2)
//      if (dot(m_Ffric,ds)> 0.0 ) {
// criterion (3)
      if ((m_Ffric+m_ks*ds).norm() > m_Ffric.norm()) { // tangential force greater than static friction -> dynamic
        //m_Ffric= m_mu_d*force.norm() *ds/ds.norm();
        m_Ffric = m_mu_d*force.norm()*(m_Ffric+m_ks*ds)/(m_Ffric+m_ks*ds).norm();
        m_force_deficit = Vec3(0.0,0.0,0.0);
        m_E_diss = fabs(m_Ffric * ds); // energy dissipated
  //      cout <<" continue slip :  " << m_Ffric<< "norm  "<< m_Ffric.norm()<< endl;
      } else { // static friction or no frictional force
        m_is_slipping= false;
        m_Ffric += m_ks*ds;
        m_E_diss = 0.0; // no energy dissipated
//        cout <<" slip-> stick:  " << m_Ffric<< endl;
      }
    }
    const Vec3 Moment1(cross(pos-m_p1->getPos(),  m_Ffric));
    const Vec3 Moment2(cross(pos-m_p2->getPos(), -m_Ffric));
    m_p1->applyForce(m_Ffric,pos);
    m_p2->applyForce(-1.0*m_Ffric,pos);
    m_p1->applyMoment(Moment1);
    m_p2->applyMoment(Moment2);
    m_cpos=pos;
    m_is_touching=true;
  } else { // no contact -> all forces are 0

//    cout << m_p1->getID() << "---  " << m_p2->getID() << "departed " <<endl;
    m_Ffric=Vec3(0.0,0.0,0.0);
    m_force_deficit=Vec3(0.0,0.0,0.0);
    m_normal_force=Vec3(0.0,0.0,0.0);
    m_is_slipping=false;
    m_is_touching=false;
    m_E_diss=0.0; // no energy dissipated
  }
}

bool CRotThermFrictionInteraction::isPersistent()
{
  Vec3 D=m_p1->getPos()-m_p2->getPos();
  double dist=D*D;
  double eq_dist=m_p1->getRad()+m_p2->getRad();
  return dist<=(eq_dist*eq_dist);
}

/*!
  get the force needed to overcome friction and make the interaction slip
*/
double CRotThermFrictionInteraction::getAbsForceDeficit()const
{
  return m_force_deficit.norm();
}

/*!
 Calculate the normal force.
*/
void CRotThermFrictionInteraction::calcNormalForce()
{
  Vec3 pos;
  // calculate distance
  const Vec3 D=m_p1->getPos()-m_p2->getPos();
  double dist=D*D;
  const double eq_dist=m_p1->getRad()+m_p2->getRad();
  // check if there is contact
  if(dist<(eq_dist*eq_dist)){ // contact -> calculate forces
    //--- elastic force ---
    dist=sqrt(dist);
    m_normal_force=D*(m_k*(dist-eq_dist)/dist);
  }
}

/*!
  get the potential energy stored in the interaction
  
  \warning For performance reasons the tangential part of the elastic energy is calculated directly from the current tangential force whereas it would be more accurate to calculate it incrementally during the force calculation. Be aware that therefore the potential energy is an approximation. Tests suggest an accuracy of a few percent.
*/
double CRotThermFrictionInteraction::getPotentialEnergy() const
{
  double e_pot_norm=0.5*m_normal_force*m_normal_force/m_k;
  double e_pot_tan=0.5*m_Ffric*m_Ffric/m_ks;

  return e_pot_norm+e_pot_tan;
}

/*!
  Get the static/dynamic status of the interaction. Returns 1 for a contact in dynamic 
  friction, 0 for static or no contact
*/
double CRotThermFrictionInteraction::getSlipping() const
{
  double res=m_is_slipping ? 1.0 : 0.0;
  return res;
}

/*!
  Get the contact status of the interaction. Returns 1 for an actual contact, 
  0 for no contact
*/
double CRotThermFrictionInteraction::getTouching() const
{
  double res=m_is_touching ? 1.0 : 0.0;
  return res;
}

/*!
  Get "sticking" contacts, i.e. return 1 if the contact is touching but not 
  slipping, 0 otherwise
*/
double CRotThermFrictionInteraction::getSticking() const
{
  const double res=(m_is_touching && !m_is_slipping) ? 1.0 : 0.0;
  return res;
}

/*!
  return the amount of energy dissipated during the last time step
*/
double CRotThermFrictionInteraction::getDissipatedEnergy() const
{
  return m_E_diss;
}

Vec3 CRotThermFrictionInteraction::getForce() const
{
  const Vec3 f=m_is_touching ? m_Ffric-m_normal_force : Vec3(0.0,0.0,0.0);
  return f; 
}

/*!
  If the particles are in contact, get normal force, if not in contact return (0,0,0)
*/
Vec3 CRotThermFrictionInteraction::getNormalForce() const
{
  const Vec3 f=m_is_touching ? m_normal_force : Vec3(0.0,0.0,0.0);
  return f; 
}

/*
calculate heat transferred bwtween 2 particles and applied to them 

*/
void CRotThermFrictionInteraction::calcHeatTrans()
{

  double  eta = 1.5 ;

  Vec3 D=m_p1->getPos()-m_p2->getPos();
  double Rij2 = D.norm2();
  double d_temp = m_p2->getTemperature() -  m_p1->getTemperature() ;
  double heatij = eta*m_diffusivity*d_temp/Rij2 ;
 
   m_p1->applyHeatTrans(heatij) ;
   m_p2->applyHeatTrans(-heatij) ;
}

/*
 calculate heat generated by driction at one time step and distributed it to two particles 
*/

void CRotThermFrictionInteraction::calcHeatFrict()
{
  double heat_frict = 0.0 ;
  double ratio,heati,heatj ;

   if(getSlipping()) {

     heat_frict = dot(m_Ffric, m_ds) ; 

     if(m_p1->getDo2dCalculations()){
        ratio = m_p1->getRad()*m_p1->getRad()/(m_p1->getRad()*m_p1->getRad()+m_p2->getRad()*m_p2->getRad() );
     }else {
        ratio = m_p1->getRad()*m_p1->getRad()*m_p1->getRad()/(m_p1->getRad()*m_p1->getRad()*m_p1->getRad()+m_p2->getRad()*m_p2->getRad()*m_p2->getRad() );
     }

     heati =  heat_frict*ratio ;
     heatj =  heat_frict*(1.0-ratio) ;
     m_p1->applyHeatFrict(heati) ;
     m_p2->applyHeatFrict(heatj) ;
   }
}


/*!
  Get the particle member function which returns a scalar field of a given name.

  \param name the name of the field 
*/
CRotThermFrictionInteraction::ScalarFieldFunction CRotThermFrictionInteraction::getScalarFieldFunction(const string& name)
{
  CRotThermFrictionInteraction::ScalarFieldFunction sf;

  if(name=="force_deficit"){
    sf=&CRotThermFrictionInteraction::getAbsForceDeficit;
  } else if (name=="potential_energy"){
    sf=&CRotThermFrictionInteraction::getPotentialEnergy;
  } else if (name=="slipping"){
    sf=&CRotThermFrictionInteraction::getSlipping;
  } else if (name=="sticking"){
    sf=&CRotThermFrictionInteraction::getSticking;
  } else if (name=="count"){
    sf=&CRotThermFrictionInteraction::Count;
  } else if (name=="dissipated_energy") {
    sf=&CRotThermFrictionInteraction::getDissipatedEnergy;
  } else {
    sf=NULL;
    cerr << "ERROR - invalid name for interaction scalar  access function" << endl; 
  }
  
  return sf;
}

CRotThermFrictionInteraction::CheckedScalarFieldFunction CRotThermFrictionInteraction::getCheckedScalarFieldFunction(const string& name)
{
  CRotThermFrictionInteraction::CheckedScalarFieldFunction sf = NULL;
  cerr << "ERROR - invalid name for interaction scalar  access function" << endl; 
  return sf;
}

/*!
  Get the particle member function which returns a vector field of a given name.

  \param name the name of the field 
*/
CRotThermFrictionInteraction::VectorFieldFunction CRotThermFrictionInteraction::getVectorFieldFunction(const string& name)
{
  CRotThermFrictionInteraction::VectorFieldFunction vf=NULL;

  if (name=="force") {
    vf = &CRotThermFrictionInteraction::getForce;    
  } else if (name=="normal_force") {
    vf = &CRotThermFrictionInteraction::getNormalForce;    
  } else {
    cerr << "ERROR - invalid name for interaction vector access function" << endl;
  }

  return vf;
}

/*!
  Pack a CFrictionInteraction into a TML packed message

  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::pack<CRotThermFrictionInteraction>(const CRotThermFrictionInteraction& I)
{
  append(I.m_k);
  append(I.m_r0);
  append(I.m_mu_d);
  append(I.m_mu_s);
  append(I.m_ks);
  append(I.m_dt);
  append(I.m_id[0]);
  append(I.m_id[1]);
  append(static_cast<int>(I.m_is_slipping));
  append(I.m_Ffric.X());
  append(I.m_Ffric.Y());
  append(I.m_Ffric.Z());

  append(I.m_diffusivity);
  append(I.m_ds.X()) ;
  append(I.m_ds.Y()) ;
  append(I.m_ds.Z()) ;
}

/*!
  Unpack a CFrictionInteraction from a TML packed message

  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::unpack<CRotThermFrictionInteraction>(CRotThermFrictionInteraction& I)
{
  I.m_k=pop_double();
  I.m_r0=pop_double();
  I.m_mu_d=pop_double();
  I.m_mu_s=pop_double();
  I.m_ks=pop_double();
  I.m_dt=pop_double();
  I.m_id.clear();
  I.m_id.push_back(pop_int());
  I.m_id.push_back(pop_int());
  I.m_is_slipping = static_cast<bool>(pop_int());
  I.m_Ffric.X() = pop_double();
  I.m_Ffric.Y() = pop_double();
  I.m_Ffric.Z() = pop_double();
  
  I.m_diffusivity = pop_double();
  I.m_ds.X() = pop_double();
  I.m_ds.Y() = pop_double();
  I.m_ds.Z() = pop_double();
}

/*!
  Save restart data to an open ostream 

  \param oStream the output stream
*/
void CRotThermFrictionInteraction::saveRestartData(std::ostream &oStream)
{
  oStream << m_k << " ";
  oStream << m_r0 << " ";
  oStream << m_mu_d << " ";
  oStream << m_mu_s << " ";
  oStream << m_ks << " ";
  oStream << m_dt << " ";
  oStream << m_id[0] << " ";
  oStream << m_id[1] << " ";
  oStream << m_is_slipping << " ";
  oStream << m_is_touching << " ";
  oStream << m_Ffric.X() << " ";
  oStream << m_Ffric.Y() << " ";
  oStream << m_Ffric.Z() << " ";
  oStream << m_diffusivity << " ";
  oStream << m_ds.X() << " ";
  oStream << m_ds.Y() << " ";
  oStream << m_ds.Z();
}


/*!
  Load restart data from an open istream 

  \param iStream the input stream
*/
void CRotThermFrictionInteraction::loadRestartData(std::istream &iStream)
{
  iStream >> m_k ;
  iStream >> m_r0 ;
  iStream >> m_mu_d ;
  iStream >> m_mu_s ;
  iStream >> m_ks ;
  iStream >> m_dt ;
  iStream >> m_id[0] ;
  iStream >> m_id[1] ;
  iStream >> m_is_slipping ;
  iStream >> m_is_touching ;
  iStream >> m_Ffric.X() ;
  iStream >> m_Ffric.Y() ;
  iStream >> m_Ffric.Z();
  iStream >> m_diffusivity;
  iStream >> m_ds.X();
  iStream >> m_ds.Y();
  iStream >> m_ds.Z();
}

ostream& operator<<(ostream& ost,const CRotThermFrictionInteraction& BI)
{
  ost << "[" << BI.m_p1->getID() << " - " << BI.m_p2->getID() << "]";
  return ost;
}
