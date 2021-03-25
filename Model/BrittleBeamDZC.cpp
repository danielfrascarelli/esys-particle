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

#include "Model/BrittleBeamDZC.h"

// --- TML includes ---
#include "tml/message/packed_message_interface.h"

#include <stdexcept>


/*!
    Default constructor for rotational bonded interaction parameters using
    stress failure criterion.  Sets all parameters to 0.0.
*/
BrittleBeamDZCIGP::BrittleBeamDZCIGP() : ARotBondedIGP(), cohesion(0.0), tCutoff(0.0), fAngle(0.0), cCutoff(0.0), beta1(0.0), beta2(0.0)
{}
    
/*!
    Constructor for rotational bonded interaction parameters using
    stress failure criterion. 
*/
BrittleBeamDZCIGP::BrittleBeamDZCIGP(
        const  std::string &name,
        double youngsModulus,
        double poissonsRatio,
        double cohesion,
        double frictionAngle,
        double tensionCutoff,
        double compressCutoff,
        double beta1,
        double beta2,
        int tag) : ARotBondedIGP(name, youngsModulus, poissonsRatio,tag),cohesion(cohesion), tCutoff(tensionCutoff), fAngle(frictionAngle), cCutoff(compressCutoff), beta1(beta1), beta2(beta2)
{}

/*!
    Constructor for rotational bonded interaction parameters using
    stress failure criterion. Using raw bond stiffnesses as input.
*/
BrittleBeamDZCIGP::BrittleBeamDZCIGP(const std::string &name,
        double kr,
        double ks,
        double kt,
        double kb,
        double cohesion,
        double frictionAngle,
        double tensionCutoff,
        double compressCutoff,
        double beta1,
        double beta2,
        int tag
    ) : ARotBondedIGP(name, kr, ks, kt, kb, tag, true) ,cohesion(cohesion), tCutoff(tensionCutoff), fAngle(frictionAngle), cCutoff(compressCutoff), beta1(beta1), beta2(beta2)
{}


// --------- END PARAMETER CLASS ----------


/*!
    Default constructor for rotational bonded interactions.
    Just sets all parameters to 0.0 
*/
BrittleBeamDZCInteraction::BrittleBeamDZCInteraction() : ARotBondedInteraction(),
    m_tCutoff(0.0),
    m_cCutoff(0.0),
    m_cohesion(0.0),
    m_tanAngle(0.0),
    m_beta1(0.0),
    m_beta2(0.0), 
    m_effR(0.0), 
    m_effA(0.0), 
    m_effI(0.0), 
    m_effJ(0.0)
{}


/*!
    Constructor for base of rotational bonded interactions using two particle pointers
    and a parameter set contained in a ARotBondedIGP.
*/
BrittleBeamDZCInteraction::BrittleBeamDZCInteraction(CRotParticle* p1,CRotParticle* p2,const BrittleBeamDZCIGP& param) : 
    ARotBondedInteraction(p1,p2,ARotBondedIGP(param)),
    m_tCutoff(param.tCutoff),
    m_cohesion(param.cohesion),
    m_tanAngle(std::tan(param.fAngle/57.2957795)),
    m_cCutoff(param.cCutoff),
    m_beta1(param.beta1),
    m_beta2(param.beta2)
{
   m_effR = 0.5*m_dist;
   m_effA = M_PI*m_effR*m_effR;
   m_effI = 0.25*M_PI*m_effR*m_effR*m_effR*m_effR;
   m_effJ = 0.5*M_PI*m_effR*m_effR*m_effR*m_effR;

}

/*!
    Default destructor - does nothing
*/
BrittleBeamDZCInteraction::~BrittleBeamDZCInteraction()
{
}

/*!
    Check if the fracture criterion has been exceeded. If so, flag the particles 
    (for the rebuilding of the other interactions) and return "true", so the update
    of this interaction group can remove the interaction.

    Applies the modified Mohr-Coulomb failure criterion of Ding and Zhang (2014)
    Criterion has compressive and tensile strength cutoffs, as well as 
    suppression of contributions of bending and torsion moments.
*/
bool BrittleBeamDZCInteraction::broken()
{
    bool res;
    double shearStress, normalStress;

    shearStress = m_shForce/m_effA + m_beta2*m_tMoment*m_effR/m_effJ;
    normalStress = -1.*(m_nForce/m_effA  + m_beta1*m_bMoment*m_effR/m_effI);
  
    // check failure criterion
    res=(normalStress < -1.*m_tCutoff) || (shearStress > m_cCutoff) || (shearStress > (m_cohesion + normalStress*m_tanAngle));
    
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
BrittleBeamDZCInteraction::ScalarFieldFunction BrittleBeamDZCInteraction::getScalarFieldFunction(const string& name)
{
  BrittleBeamDZCInteraction::ScalarFieldFunction sf;
                                                                                
  if (name=="potential_energy"){
    sf=&BrittleBeamDZCInteraction::getPotentialEnergy;
  } else if (name=="e_pot_normal"){
    sf=&BrittleBeamDZCInteraction::getNormalPotentialEnergy;
  } else if (name=="e_pot_shear"){
    sf=&BrittleBeamDZCInteraction::getShearPotentialEnergy;
  } else if (name=="e_pot_twist"){
    sf=&BrittleBeamDZCInteraction::getTwistPotentialEnergy;
  } else if (name=="e_pot_bend"){
    sf=&BrittleBeamDZCInteraction::getBendPotentialEnergy;
  } else if (name=="count"){
    sf=&BrittleBeamDZCInteraction::Count;
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
BrittleBeamDZCInteraction::VectorFieldFunction BrittleBeamDZCInteraction::getVectorFieldFunction(const string& name)
{
  BrittleBeamDZCInteraction::VectorFieldFunction vf;
                                                                                
  if (name=="force"){
    vf=&BrittleBeamDZCInteraction::getForce;
  } else if (name=="normal_force"){
    vf=&BrittleBeamDZCInteraction::getNormalForce;
  } else if (name=="tangential_force"){
    vf=&BrittleBeamDZCInteraction::getTangentialForce;
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
BrittleBeamDZCInteraction::CheckedScalarFieldFunction BrittleBeamDZCInteraction::getCheckedScalarFieldFunction(const string& name)
{
  BrittleBeamDZCInteraction::CheckedScalarFieldFunction sf;

  sf=NULL;
  cerr << "ERROR - invalid name for interaction scalar  access function" << endl;

  return sf;
}

void BrittleBeamDZCInteraction::saveCheckPointData(std::ostream &oStream)
{
  throw std::runtime_error("BrittleBeamDZCInteraction::saveCheckPointData not implemented.");
}

void BrittleBeamDZCInteraction::loadCheckPointData(std::istream &iStream)
{
  throw std::runtime_error("BrittleBeamDZCInteraction::loadCheckPointData not implemented.");
}

/*!
  Pack this object into a TML packed message

  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::pack<BrittleBeamDZCInteraction>(const BrittleBeamDZCInteraction& I)
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
    append(I.m_tCutoff);
    append(I.m_cCutoff);
    append(I.m_cohesion);
    append(I.m_tanAngle);
    append(I.m_beta1);
    append(I.m_beta2);
}

/*!
  Unpack a CBondedInteraction from a TML packed message
                                                                                
  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::unpack<BrittleBeamDZCInteraction>(BrittleBeamDZCInteraction& I)
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
    I.m_tCutoff = pop_double();
    I.m_cCutoff = pop_double();
    I.m_cohesion = pop_double();
    I.m_tanAngle = pop_double();
    I.m_beta1 = pop_double();
    I.m_beta2 = pop_double();
}

/*!
  save restart data to ostream

  \param oStream the output stream
*/
void BrittleBeamDZCInteraction::saveRestartData(std::ostream &oStream)
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
    oStream << m_tCutoff << " ";
    oStream << m_cCutoff << " ";
    oStream << m_cohesion << " ";
    oStream << m_tanAngle << " ";
    oStream << m_beta1 << " ";
    oStream << m_beta2 << std::endl;
}

/*!
  load restart data from stream

  \param iStream the input stream
*/
void BrittleBeamDZCInteraction::loadRestartData(std::istream &iStream)
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
    iStream >> m_tCutoff;
    iStream >> m_cCutoff;
    iStream >> m_cohesion;
    iStream >> m_tanAngle;
    iStream >> m_beta1;
    iStream >> m_beta2;
    
    setTag(tag);
}
