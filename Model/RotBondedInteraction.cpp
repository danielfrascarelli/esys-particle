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

#include "Model/RotBondedInteraction.h"
#include "Model/BondedInteractionCpData.h"
#include "Foundation/console.h"

// --- TML includes ---
#include "tml/message/packed_message_interface.h"

#include <stdexcept>


/*!
    Default constructor for basic rotational bonded interaction parameters.
    Sets all stiffnesses to 0.0
*/
CRotBondedIGP::CRotBondedIGP() : ARotBondedIGP(),
   max_nForce(0.0),
   max_shForce(0.0),
   max_tMoment(0.0),
   max_bMoment(0.0),
   meanR_scaling(true),
   truncated(1.0),
   beta1(1.0),
   beta2(1.0)
{}

/*!
    Constructor for basic rotational bonded interaction parameters from
    individual stiffness parameters (Wang et al 2006)
*/
CRotBondedIGP::CRotBondedIGP(
  const  std::string &name,
  double kr,
  double ks,
  double kt,
  double kb,
  double max_nForce,
  double max_shForce,
  double max_tMoment,
  double max_bMoment,
  int    tag,
  bool   scaling,
  bool   AmeanR_scaling,
  double truncated,
  double beta1,
  double beta2
) : ARotBondedIGP(name,kr,ks,kt,kb,tag,scaling),
   max_nForce(max_nForce),
   max_shForce(max_shForce),
   max_tMoment(max_tMoment),
   max_bMoment(max_bMoment),
   meanR_scaling(AmeanR_scaling),
   truncated(truncated),
   beta1(beta1),
   beta2(beta2)
{}

/*!
    Constructor for basic rotational bonded interaction parameters from
    brittle beam parameters (stiffness only, not strength).
*/
CRotBondedIGP::CRotBondedIGP(
  const  std::string &name,
  double youngsModulus,
  double poissonsRatio,
  double cohesion,
  double tanAngle,
  int    tag,
  bool AmeanR_scaling,
  double truncated,
  double beta1,
  double beta2
) : ARotBondedIGP(name,youngsModulus,poissonsRatio,tag),
   meanR_scaling(AmeanR_scaling),
   truncated(truncated),
   beta1(beta1),
   beta2(beta2)
{
   max_nForce = M_PI*cohesion/tanAngle;
   max_shForce = M_PI*cohesion;
   max_bMoment = M_PI*cohesion/tanAngle;
   max_tMoment = M_PI*cohesion;
}


// --------- END PARAMETER CLASS ----------


/*!
    Default constructor for base of rotational bonded interactions.
    Just sets all parameters to 0.0 
*/
CRotBondedInteraction::CRotBondedInteraction():ARotBondedInteraction()
{
  m_max_nForce  = 0.0 ;
  m_max_shForce = 0.0 ;
  m_max_tMoment = 0.0 ;
  m_max_bMoment = 0.0 ;

  m_meanR_scaling = true;
  m_truncated = 1.0;
  m_beta1 = 1.0;
  m_beta2 = 1.0;
}


/*!
    Constructor for base of rotational bonded interactions using two particle pointers
    and a parameter set contained in a ARotBondedIGP.
*/
CRotBondedInteraction::CRotBondedInteraction(CRotParticle* p1,CRotParticle* p2,const CRotBondedIGP& param):ARotBondedInteraction(p1,p2,ARotBondedIGP(param))
{
  double effR=1.0;
  double effL=1.0;
  double effA=1.0;
  double momentJ=1.0;
  double momentI=1.0;
  m_scaling = param.scaling;
  m_meanR_scaling = param.meanR_scaling;
  // scale elastic param
  if (m_scaling == true) {
    if(!CParticle::getDo2dCalculations()){

      if (m_meanR_scaling == true) {
        effR=0.5*m_dist;
      }
      else {
        effR=fmin(p1->getRad(), p2->getRad());
      }
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
  
  m_max_nForce  = param.max_nForce * effA;
  m_max_shForce = param.max_shForce * effA;
  m_max_bMoment = param.max_bMoment * momentI / effR;
  m_max_tMoment = param.max_tMoment * momentJ / effR; 

  m_truncated = param.truncated;
  m_beta1 = param.beta1;
  m_beta2 = param.beta2;
}

/*!
    Default destructor - does nothing
*/
CRotBondedInteraction::~CRotBondedInteraction()
{
}

/*!
  Check if the fracture criterion has been exceeded. If so, flag the particles 
  (for the rebuilding of the other interactions) and return "true", so the update
  of this interaction group can remove the interaction.
*/
bool CRotBondedInteraction::broken()
{
  double shearStress, normalStress;
  bool res;
/*
// Original failure criterion:
  const double crit_nf = (m_nForce/m_max_nForce > 0.0) ? m_nForce/m_max_nForce :0.0;
  const double criterion =
    crit_nf  +
    m_shForce/m_max_shForce +
    m_tMoment/m_max_tMoment +
    m_bMoment/m_max_bMoment;
*/

/*
// Mohr-Coulomb-like failure criterion:
  const double criterion =
    m_nForce/m_max_nForce  +
    m_shForce/m_max_shForce +
    m_tMoment/m_max_tMoment +
    m_bMoment/m_max_bMoment;


  if(criterion > 1.0) { // broken
    if (m_p1!=NULL) m_p1->setFlag();
    if (m_p2!=NULL) m_p2->setFlag();
    res=true;
  } else {
    res=false;
  }
*/

// direct Mohr-Coulomb failure criterion:
  shearStress = m_shForce/m_max_shForce + m_beta2*m_tMoment/m_max_tMoment ;
  normalStress = m_nForce/m_max_nForce  + m_beta1*m_bMoment/m_max_bMoment ;

  if ((normalStress > m_truncated) || (shearStress > (1.0 - normalStress))) { // broken
    if (m_p1!=NULL) m_p1->setFlag();
    if (m_p2!=NULL) m_p2->setFlag();
    res=true;
  } else {
    res=false;
  }
  return res;
}

double CRotBondedInteraction::getCriterion() const
{
  return m_nForce/m_max_nForce   +
    m_shForce/m_max_shForce +
    m_tMoment/m_max_tMoment +
    m_bMoment/m_max_bMoment;
}


/*!
  Get the particle member function which returns a scalar field of a given name.

  \param name the name of the field 
*/
CRotBondedInteraction::ScalarFieldFunction CRotBondedInteraction::getScalarFieldFunction(const string& name)
{
  CRotBondedInteraction::ScalarFieldFunction sf;
                                                                                
  if (name=="potential_energy"){
    sf=&CRotBondedInteraction::getPotentialEnergy;
  } else if (name=="e_pot_normal"){
    sf=&CRotBondedInteraction::getNormalPotentialEnergy;
  } else if (name=="e_pot_shear"){
    sf=&CRotBondedInteraction::getShearPotentialEnergy;
  } else if (name=="e_pot_twist"){
    sf=&CRotBondedInteraction::getTwistPotentialEnergy;
  } else if (name=="e_pot_bend"){
    sf=&CRotBondedInteraction::getBendPotentialEnergy;
  } else if (name=="count"){
    sf=&CRotBondedInteraction::Count;
  } else if (name=="breaking_criterion"){
    sf=&CRotBondedInteraction::getCriterion;
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
CRotBondedInteraction::VectorFieldFunction CRotBondedInteraction::getVectorFieldFunction(const string& name)
{
  CRotBondedInteraction::VectorFieldFunction vf;
                                                                                
  if (name=="force"){
    vf=&CRotBondedInteraction::getForce;
  } else if (name=="normal_force"){
    vf=&CRotBondedInteraction::getNormalForce;
  } else if (name=="tangential_force"){
    vf=&CRotBondedInteraction::getTangentialForce;
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
CRotBondedInteraction::CheckedScalarFieldFunction CRotBondedInteraction::getCheckedScalarFieldFunction(const string& name)
{
  CRotBondedInteraction::CheckedScalarFieldFunction sf;

  sf=NULL;
  cerr << "ERROR - invalid name for interaction scalar  access function" << endl;

  return sf;
}

void CRotBondedInteraction::saveCheckPointData(std::ostream &oStream)
{
  BondedInteractionCpData(*this).saveCheckPointData(oStream);
}

void CRotBondedInteraction::loadCheckPointData(std::istream &iStream)
{
  throw std::runtime_error("CRotBondedInteraction::loadCheckPointData not implemented.");
}

/*!
  Pack this object into a TML packed message

  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::pack<CRotBondedInteraction>(const CRotBondedInteraction& I)
{
  append(I.m_id[0]);
  append(I.m_id[1]);
  append(static_cast<int>(I.m_init));
  append(I.m_kr);
  append(I.m_ks);
  append(I.m_kb);
  append(I.m_kt);
  append(static_cast<int>(I.m_scaling));
  append(static_cast<int>(I.m_meanR_scaling));
  append(I.m_truncated);
  append(I.m_beta1);
  append(I.m_beta2);
  append(I.m_nForce);  // exchange of current forces is needed 
  append(I.m_shForce); // to ensure the bond is broken correctly
  append(I.m_tMoment); // in case of a calc - exchange - update (i.e. break)
  append(I.m_bMoment); // sequence because broken() needs the forces
  append(I.m_max_nForce);
  append(I.m_max_shForce);
  append(I.m_max_tMoment);
  append(I.m_max_bMoment);
  append(I.m_D);
  append(I.m_tag);
}

/*!
  Unpack a CBondedInteraction from a TML packed message
                                                                                
  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::unpack<CRotBondedInteraction>(CRotBondedInteraction& I)
{
  I.m_id.clear();
  I.m_id.push_back(pop_int());
  I.m_id.push_back(pop_int());
  I.m_init = static_cast<bool>(pop_int());
  I.m_kr = pop_double();
  I.m_ks = pop_double();
  I.m_kb = pop_double();
  I.m_kt = pop_double();
  I.m_scaling=static_cast<bool>(pop_int());
  I.m_meanR_scaling=static_cast<bool>(pop_int());
  I.m_truncated = pop_double();
  I.m_beta1 = pop_double();
  I.m_beta2 = pop_double();
  I.m_nForce = pop_double();
  I.m_shForce = pop_double();
  I.m_tMoment = pop_double();
  I.m_bMoment = pop_double();
  I.m_max_nForce = pop_double();
  I.m_max_shForce = pop_double();
  I.m_max_tMoment = pop_double();
  I.m_max_bMoment = pop_double();
  I.m_D=pop_vec3();
  I.m_tag=pop_int();
}

/*!
  save restart data to ostream

  \param oStream the output stream
*/
void CRotBondedInteraction::saveRestartData(std::ostream &oStream)
{
  oStream << m_id[0] << " ";
  oStream << m_id[1] << " ";
  oStream << m_init  << " ";
  oStream << m_kr << " ";
  oStream << m_ks << " ";
  oStream << m_kb << " ";
  oStream << m_kt << " ";
  oStream << m_scaling << " ";
  oStream << m_meanR_scaling << " ";
  oStream << m_truncated << " ";
  oStream << m_beta1 << " ";
  oStream << m_beta2 << " ";
  oStream << m_max_nForce << " ";
  oStream << m_max_shForce << " ";
  oStream << m_max_tMoment << " ";
  oStream << m_max_bMoment << " ";
  oStream << m_D << " ";
  oStream << getTag() << std::endl;
}

/*!
  load restart data from stream

  \param iStream the input stream
*/
void CRotBondedInteraction::loadRestartData(std::istream &iStream)
{
  int tag;
  iStream >> m_id[0];
  iStream >> m_id[1];
  iStream >> m_init;
  iStream >> m_kr;
  iStream >> m_ks;
  iStream >> m_kb;
  iStream >> m_kt;
  iStream >> m_scaling;
  iStream >> m_meanR_scaling;
  iStream >> m_truncated;
  iStream >> m_beta1;
  iStream >> m_beta2;
  iStream >> m_max_nForce;
  iStream >> m_max_shForce;
  iStream >> m_max_tMoment;
  iStream >> m_max_bMoment;
  iStream >> m_D;
  iStream >> tag;
  setTag(tag);
}

ostream& operator<<(ostream& ost,const CRotBondedInteraction& BI)
{
  ost << "[" << BI.m_p1->getID() << " - " << BI.m_p2->getID() << "]";
  return ost;
}
