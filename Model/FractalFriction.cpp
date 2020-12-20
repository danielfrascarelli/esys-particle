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

#include "Model/FractalFriction.h"
#include "Foundation/console.h"
#include "tml/message/packed_message_interface.h"

FractalFrictionIGP::FractalFrictionIGP()
  : AIGParam(),
   k(0.0),
   mu_0(0.0),
   k_s(0.0),
   dt(0.0),
   mu(),          //!< pointer to the array of friction coeff.
   x0(0.0),y0(0.0),dx(0.0),dy(0.0),  //!< origin and grid spacing of the array
   nx(0),ny(0)           //!< array size
{
}

FractalFrictionIGP::FractalFrictionIGP(const FractalFrictionIGP& F) : AIGParam(F)
{
  k    = F.k;
  mu_0 = F.mu_0;
  k_s  = F.k_s;
  dt   = F.dt;
  x0   = F.x0;
  y0   = F.y0;
  dx   = F.dx;
  dy   = F.dy;
  nx   = F.nx;
  ny   = F.ny;
  mu   = boost::shared_ptr<double>(new double[nx*ny]);
  
  for (int i = 0; i < nx*ny; i++)
  {
    (mu.get())[i] = (F.mu.get())[i];
  }
}

FractalFrictionIGP &FractalFrictionIGP::operator=(const FractalFrictionIGP &F)
{
  k    = F.k;
  mu_0 = F.mu_0;
  k_s  = F.k_s;
  dt   = F.dt;
  x0   = F.x0;
  y0   = F.y0;
  dx   = F.dx;
  dy   = F.dy;
  nx   = F.nx;
  ny   = F.ny;
  mu   = boost::shared_ptr<double>(new double[nx*ny]);
  
  for (int i = 0; i < nx*ny; i++)
  {
    (mu.get())[i] = (F.mu.get())[i];
  }

  return *this;
}

FractalFrictionIGP::~FractalFrictionIGP()
{
  ///////delete mu;
}


CFractalFriction::CFractalFriction()
{
  m_k=0.0;
  m_mu=0.0;
  m_r0=0.0;
  m_ks=0.0;
  m_dt=0.0;
}

/*!
  Construct a CFractalFriction from 2 particle pointers and the parameters

  \param p1 pointer to the first particle
  \param p2 pointer to the second particle
  \param param the interaction parameters
*/
CFractalFriction::CFractalFriction(CParticle* p1,CParticle* p2,const FractalFrictionIGP& param):CFrictionInteraction(p1,p2)
{
  m_k=param.k;
  m_ks=param.k_s;
  m_r0=p1->getRad()+p2->getRad();
  m_dt=param.dt;
  m_cpos=p1->getPos()+((p2->getPos()-p1->getPos())*p1->getRad()/m_r0);
  int idx=int(floor((m_cpos.X()-param.x0)/param.dx)); // where in the mu-array are we (x-dir)
  idx= (idx<0) ? 0 : idx;                           // clamp idx to 0..nx-1
  idx= (idx>=param.nx) ? param.nx-1 : idx;
  int idy=int(floor((m_cpos.Y()-param.y0)/param.dy)); // where in the mu-array are we (x-dir)
  idy= (idy<0) ? 0 : idy;                           // clamp idy to 0..ny-1
  idy= (idy>=param.ny) ? param.ny-1 : idy;
  int index=idx*param.ny+idy;
  m_mu=param.mu_0*param.mu.get()[index];
  //  console.Debug() << "cpos,ix,iy,index,mu: " << m_cpos << " , " << idx << " , " << idy << " , " << index << " , " << m_mu << "\n";
}

CFractalFriction::~CFractalFriction()
{}

/*!
  Get the particle member function which returns a scalar field of a given name.

  \param name the name of the field 
*/
CFractalFriction::ScalarFieldFunction CFractalFriction::getScalarFieldFunction(const string& name)
{
  CFractalFriction::ScalarFieldFunction sf;

  if (name=="potential_energy"){
    sf=&CFractalFriction::getPotentialEnergy;
  } else if (name=="slipping"){
    sf=&CFractalFriction::getSlipping;
  } else if (name=="count"){
    sf=&CFractalFriction::Count;
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
CFractalFriction::VectorFieldFunction CFractalFriction::getVectorFieldFunction(const string& name)
{
  CFractalFriction::VectorFieldFunction vf;

  vf=NULL;
  cerr << "ERROR - invalid name for interaction vector access function" << endl; 

  return vf;
}

/*!
  Get the particle member function which returns a checked scalar field of a given name.

  \param name the name of the field
*/
CFractalFriction::CheckedScalarFieldFunction CFractalFriction::getCheckedScalarFieldFunction(const string& name)
{
	CFractalFriction::CheckedScalarFieldFunction sf;

	if (name=="mu_eff_xy"){
		sf=&CFractalFriction::getMuEffXY;
	} else if (name=="mu_eff_xz"){
		sf=&CFractalFriction::getMuEffXZ;
	} else if (name=="F_fric"){
		sf=&CFractalFriction::getAbsFrictionalForce;
	} else if (name=="muF_n") {
		sf=&CFractalFriction::getAbsMuFN;
	} else if (name=="v_slip") {
		sf=&CFractalFriction::getSlipVelocity;
	} else {
		sf=NULL;
		cerr << "ERROR - invalid name for interaction scalar  access function" << endl;
	}

	return sf;
}

/*!
  Pack a CFractalFriction into a TML packed message

  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::pack<CFractalFriction>(const CFractalFriction& I)
{
  append(I.m_k);
  append(I.m_r0);
  append(I.m_mu);
  append(I.m_ks);
  append(I.m_dt);
  append(I.m_id[0]);
  append(I.m_id[1]);
}

/*!
  Unpack a CFractalFriction from a TML packed message

  \param I the interaction
*/
template<>
void TML_PackedMessageInterface::unpack<CFractalFriction>(CFractalFriction& I)
{
  I.m_k=pop_double();
  I.m_r0=pop_double();
  I.m_mu=pop_double();
  I.m_ks=pop_double();
  I.m_dt=pop_double();
  I.m_id.erase(I.m_id.begin(),I.m_id.end());
  I.m_id.push_back(pop_int());
  I.m_id.push_back(pop_int());
}

ostream& operator<<(ostream& ost,const CFractalFriction& FI)
{
  ost << "[" << FI.m_p1->getID() << " - " << FI.m_p2->getID() << "] : " << FI.m_mu ;
  return ost;
}
