/////////////////////////////////////////////////////////////
//                                                         //
// Copyright (c) 2003-2019 by The University of Queensland //
// Centre for Geoscience Computing                         //
// http://earth.uq.edu.au/centre-geoscience-computing      //
//                                                         //
// Primary Business: Brisbane, Queensland, Australia       //
// Licensed under the Open Software License version 3.0    //
// http://www.apache.org/licenses/LICENSE-2.0              //
//                                                         //
/////////////////////////////////////////////////////////////

#include "Model/BrittleBeamSC.h"

// --- TML includes ---
#include "tml/message/packed_message_interface.h"

#include <stdexcept>


/*!
    Default constructor for rotational bonded interaction parameters using
    stress failure criterion.  Sets all parameters to 0.0.
*/
BrittleBeamSCIGP::BrittleBeamSCIGP() : ARotBondedIGP(), cohesion(0.0), tCutoff(0.0), fAngle(0.0)
{}
    
/*!
    Constructor for rotational bonded interaction parameters using
    stress failure criterion. 
*/
BrittleBeamSCIGP::BrittleBeamSCIGP(
        const  std::string &name,
        double youngsModulus,
        double poissonsRatio,
        double cohesion,
        double frictionAngle,
        double tensionCutoff,
        int tag) : ARotBondedIGP(name, youngsModulus, poissonsRatio,tag),cohesion(cohesion), tCutoff(tensionCutoff), fAngle(frictionAngle)
{}

/*!
    Constructor for rotational bonded interaction parameters using
    stress failure criterion. Using raw bond stiffnesses as input.
*/
BrittleBeamSCIGP::BrittleBeamSCIGP(const std::string &name,
        double kr,
        double ks,
        double kt,
        double kb,
        double cohesion,
        double frictionAngle,
        double tensionCutoff,
        int tag
    ) : ARotBondedIGP(name, kr, ks, kt, kb, tag, true) ,cohesion(cohesion), tCutoff(tensionCutoff), fAngle(frictionAngle)
{}


// --------- END PARAMETER CLASS ----------


/*!
    Default constructor for rotational bonded interactions.
    Just sets all parameters to 0.0 
*/
BrittleBeamSCInteraction::BrittleBeamSCInteraction() : ARotBondedInteraction(),
    m_sigma3(0.0),
    m_cohesion(0.0),
    m_tanAngle(0.0) 
{}


/*!
    Constructor for base of rotational bonded interactions using two particle pointers
    and a parameter set contained in a ARotBondedIGP.
*/
BrittleBeamSCInteraction::BrittleBeamSCInteraction(CRotParticle* p1,CRotParticle* p2,const BrittleBeamSCIGP& param) : 
    ARotBondedInteraction(p1,p2,ARotBondedIGP(param)),
    m_sigma3(param.tCutoff),
    m_cohesion(param.cohesion),
    m_tanAngle(std::tan(param.fAngle/57.2957795)) 
{}

/*!
    Default destructor - does nothing
*/
BrittleBeamSCInteraction::~BrittleBeamSCInteraction()
{
}

/*!
    Check if the fracture criterion has been exceeded. If so, flag the particles 
    (for the rebuilding of the other interactions) and return "true", so the update
    of this interaction group can remove the interaction.

    Applies a truncated Mohr-Coulomb criterion to the rotated average particle stress, i.e.
    - average particle stresses
    - extract symmetric part
    - get eigenvalues 
    - get sigma_1 & sigma_3 from min/max ev
    - apply M-C in sigma_1 - sigma_3 space
*/
bool BrittleBeamSCInteraction::broken()
{
    bool res;
    Vec3 n,s,t;
    
    // get bond orientation from particle position -> use as z-axis of local coordiante system
    n=(m_p2->getPos()-m_p1->getPos()).unit();
        
    // get 2nd axis
    if (n!=Vec3(1.0,0.0,0.0)){
        s=cross(n,Vec3(1.0,0.0,0.0)).unit();
    } else {
        s=cross(n,Vec3(0.0,1.0,0.0)).unit();
    }
    
    // 3rd axis
    t=cross(n,s);
    
    // transformation matrix
    Matrix3 tmat(n,s,t);
    tmat=tmat.trans();
     
    // get stresses from the two particles, calculate average & get symmetric part
    Matrix3 sigma_p1=m_p1->sigma();
    Matrix3 sigma_p2=m_p2->sigma();
    Matrix3 sigma_avg=0.5*(sigma_p1+sigma_p2);
    Matrix3 sigma_sym=sigma_avg.get_symmetric_part();

    // rotate stress to bond-local coordinate system 
    Matrix3 sigma_rot=tmat.dot(sigma_sym).dot(tmat.trans());
   
    // extract normal and shear stress 
    double sigma_n=sigma_rot(0,0);
    double tau_x=sigma_rot(0,1);
    double tau_y=sigma_rot(0,2);
    double tau=std::sqrt(tau_x*tau_x+tau_y*tau_y);
    
    // check failure criterion - s3 < tensile cutoff || s1 > MC
    res=(sigma_n < -1.0*m_sigma3) || (tau > m_cohesion+sigma_n*m_tanAngle);
    if (res) {
        if (m_p1!=NULL) m_p1->setFlag();
        if (m_p2!=NULL) m_p2->setFlag();
    }
    
    return res;
}

/*!
  Get the particle member function which returns a scalar field of a given name.

  \param name the name of the field 
*/
BrittleBeamSCInteraction::ScalarFieldFunction BrittleBeamSCInteraction::getScalarFieldFunction(const string& name)
{
  BrittleBeamSCInteraction::ScalarFieldFunction sf;
                                                                                
  if (name=="potential_energy"){
    sf=&BrittleBeamSCInteraction::getPotentialEnergy;
  } else if (name=="e_pot_normal"){
    sf=&BrittleBeamSCInteraction::getNormalPotentialEnergy;
  } else if (name=="e_pot_shear"){
    sf=&BrittleBeamSCInteraction::getShearPotentialEnergy;
  } else if (name=="e_pot_twist"){
    sf=&BrittleBeamSCInteraction::getTwistPotentialEnergy;
  } else if (name=="e_pot_bend"){
    sf=&BrittleBeamSCInteraction::getBendPotentialEnergy;
  } else if (name=="count"){
    sf=&BrittleBeamSCInteraction::Count;
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
BrittleBeamSCInteraction::VectorFieldFunction BrittleBeamSCInteraction::getVectorFieldFunction(const string& name)
{
  BrittleBeamSCInteraction::VectorFieldFunction vf;
                                                                                
  if (name=="force"){
    vf=&BrittleBeamSCInteraction::getForce;
  } else if (name=="normal_force"){
    vf=&BrittleBeamSCInteraction::getNormalForce;
  } else if (name=="tangential_force"){
    vf=&BrittleBeamSCInteraction::getTangentialForce;
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
BrittleBeamSCInteraction::CheckedScalarFieldFunction BrittleBeamSCInteraction::getCheckedScalarFieldFunction(const string& name)
{
  BrittleBeamSCInteraction::CheckedScalarFieldFunction sf;

  sf=NULL;
  cerr << "ERROR - invalid name for interaction scalar  access function" << endl;

  return sf;
}

void BrittleBeamSCInteraction::saveCheckPointData(std::ostream &oStream)
{
  throw std::runtime_error("BrittleBeamSCInteraction::saveCheckPointData not implemented.");
}

void BrittleBeamSCInteraction::loadCheckPointData(std::istream &iStream)
{
  throw std::runtime_error("BrittleBeamSCInteraction::loadCheckPointData not implemented.");
}

/*!
  Pack this object into a TML packed message

  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::pack<BrittleBeamSCInteraction>(const BrittleBeamSCInteraction& I)
{
    append(I.m_id[0]);
    append(I.m_id[1]);
    append(static_cast<int>(I.m_init));
    append(I.m_kr);
    append(I.m_ks);
    append(I.m_kb);
    append(I.m_kt);
    append(static_cast<int>(I.m_scaling));
    append(I.m_nForce);  // exchange of current forces is needed 
    append(I.m_shForce); // to ensure the bond is broken correctly
    append(I.m_tMoment); // in case of a calc - exchange - update (i.e. break)
    append(I.m_bMoment); // sequence because broken() needs the forces
    append(I.m_D);
    append(I.m_tag);
    append(I.m_sigma3);
    append(I.m_cohesion);
    append(I.m_tanAngle);
}

/*!
  Unpack a CBondedInteraction from a TML packed message
                                                                                
  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::unpack<BrittleBeamSCInteraction>(BrittleBeamSCInteraction& I)
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
    I.m_nForce = pop_double();
    I.m_shForce = pop_double();
    I.m_tMoment = pop_double();
    I.m_bMoment = pop_double();
    I.m_D=pop_vec3();
    I.m_tag=pop_int();
    I.m_sigma3 = pop_double();
    I.m_cohesion = pop_double();
    I.m_tanAngle = pop_double();
}

/*!
  save restart data to ostream

  \param oStream the output stream
*/
void BrittleBeamSCInteraction::saveRestartData(std::ostream &oStream)
{
    oStream << m_id[0] << " ";
    oStream << m_id[1] << " ";
    oStream << m_init  << " ";
    oStream << m_kr << " ";
    oStream << m_ks << " ";
    oStream << m_kb << " ";
    oStream << m_kt << " ";
    oStream << m_scaling << " ";
    oStream << m_D << " ";
    oStream << getTag() << " ";
    oStream << m_sigma3 << " ";
    oStream << m_cohesion << " ";
    oStream << m_tanAngle;
}

/*!
  load restart data from stream

  \param iStream the input stream
*/
void BrittleBeamSCInteraction::loadRestartData(std::istream &iStream)
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
    iStream >> m_D;
    iStream >> tag;
    iStream >> m_sigma3;
    iStream >> m_cohesion;
    iStream >> m_tanAngle;
    
    setTag(tag);
}
