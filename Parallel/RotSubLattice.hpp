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

#include "Model/BrittleBeamSC.h"
#include "Model/BrittleBeamDZC.h"

/*!
  Construct RotSubLattice. Calls constructor of base class.
 
  \param param Lattice parameters
  \param rank the MPI rank
  \param comm the MPI communicator
*/
template <class T>
TRotSubLattice<T>::TRotSubLattice(const esys::lsm::CLatticeParam &prm, int rank, MPI_Comm comm, MPI_Comm worker_comm):
  TSubLattice<T>(prm,rank,comm,worker_comm)
{}

/*!
  Destructor
*/
template <class T>
TRotSubLattice<T>::~TRotSubLattice()
{}

/*!
  Set the angular velocity of a particle. Parameters are received from master. 
*/
template <class T>
void TRotSubLattice<T>::setParticleAngularVelocity()
{ 
  console.Debug() << "TSubLattice<T>::setParticleAngularVelocity()\n";
  CVarMPIBuffer buffer(this->m_comm);

  buffer.receiveBroadcast(0); // get data from master
  int id=buffer.pop_int();
  Vec3 mv=buffer.pop_vector();
  this->m_ppa->forParticle(id,(void (T::*)(Vec3))(&T::setAngVel),mv);
  console.XDebug() << "end TSubLattice<T>::setParticleAngularVelocity()\n";
}

/*!
  do the actual work adding the pair interaction group (PIG)

  \param name the name of the PIG
  \param type the type of the PIG
  \param param_buffer the buffer containing the rest of the parameters
*/
template <class T>
bool TRotSubLattice<T>::doAddPIG(const string& name,const string& type,CVarMPIBuffer& param_buffer,bool tagged)
{
  bool res;

  // try interactions defined in the base class
  res=TSubLattice<T>::doAddPIG(name,type,param_buffer,tagged);
  if(!res){ // if not successfull, try interactions defined here
    
  }
  return res;

}

/*!
  Do the work for adding the damping 

  \param type the type of damping 
  \param param_buffer the buffer containing the parameters
*/
template <class T>
bool TRotSubLattice<T>::doAddDamping(const string& type,CVarMPIBuffer& param_buffer)
{
  bool res=false;

  AParallelInteractionStorage* DG=NULL;
  string damping_name;
  if (type=="RotDamping") {
    CDampingIGP *params=extractDampingIGP(&param_buffer);
    console.Debug() << "add rotational damping\n";
    DG=new ParallelInteractionStorage_Single<T,CRotDamping<T> >(this->m_ppa,*params); // dodgy
    damping_name="RotDamping";
    res=true;
    // add to map
    this->m_damping.insert(make_pair(damping_name,DG));
    this->m_damping[damping_name]->update();
  }
  else if (type=="RotLocalDamping") {
    CLocalDampingIGP *params=extractLocalDampingIGP(&param_buffer);
    console.Debug() << "add rotational damping\n";
    DG=new ParallelInteractionStorage_Single<T,CRotLocalDamping<T> >(this->m_ppa,*params); // dodgy
    damping_name="RotLocalDamping";
    res=true;
    // add to map
    this->m_damping.insert(make_pair(damping_name,DG));
    this->m_damping[damping_name]->update();
  } else {
    res=TSubLattice<T>::doAddDamping(type,param_buffer);
  }

  return res;
}

/*!
  Add bonded interaction group to the lattice. Receive the parameters from master.
  The bonds are created from the neighbor table.
*/
template <class T>
void TRotSubLattice<T>::addRotBondedIG()
{
  console.XDebug()  << "TSubLattice<T>::addRotBondedIG()\n";
  CVarMPIBuffer param_buffer(this->m_comm);
  vector<int> conns;

  // get params
  param_buffer.receiveBroadcast(0);
  int tag=param_buffer.pop_int();
  string name=string(param_buffer.pop_string());
  double kr = param_buffer.pop_double();
  double ks = param_buffer.pop_double();
  double kt = param_buffer.pop_double();
  double kb = param_buffer.pop_double();
  double max_nForce  = param_buffer.pop_double();
  double max_shForce = param_buffer.pop_double();
  double max_tMoment = param_buffer.pop_double();
  double max_bMoment = param_buffer.pop_double();
  bool   scaling = static_cast<bool>(param_buffer.pop_int());
  bool   meanR_scaling = static_cast<bool>(param_buffer.pop_int());
  double truncated = param_buffer.pop_double();
  
  conns = TSubLattice<T>::m_temp_conn[tag];

  console.XDebug()
    << "Got RotBondedIG parameters: tag=" << tag
    << ", name=" << name.c_str()
    << ", kr=" << kr
    << ", ks=" << ks
    << ", kt=" << kt
    << ", kb=" << kb
    << ", nFrc=" << max_nForce
    << ", sFrc=" << max_shForce
    << ", tMom=" << max_tMoment
    << ", bMom=" << max_bMoment
    << ", scaling=" << scaling
    << "\n";
  // setup InteractionGroup
  CRotBondedIGP param;
  param.tag=tag;
  param.kr=kr;
  param.ks=ks;
  param.kt=kt;
  param.kb=kb;
  param.max_nForce = max_nForce ;
  param.max_shForce = max_shForce ;
  param.max_tMoment = max_tMoment ;
  param.max_bMoment = max_bMoment ;
  param.scaling = scaling;
  param.meanR_scaling = meanR_scaling;
  param.truncated = truncated;
  ParallelInteractionStorage_EB<CRotParticle,CRotBondedInteraction> *B_PIS = 
    new ParallelInteractionStorage_EB<CRotParticle,CRotBondedInteraction>(this->m_ppa,param);

 // recv broadcast connection data
 // console.XDebug() << "rank=" << this->m_tml_comm.rank() << "TSubLattice<T>::addRotBondedIG(): receiving conn_data.\n";
 // vector<int> conn_data;
 // this->m_tml_comm.recv_broadcast_cont(conn_data,0);
 // console.XDebug() << "rank=" <<this->m_tml_comm.rank() << "TSubLattice<T>::addRotBondedIG(): conn_data.size()=" << conn_data.size() << "\n";
 vector<int> vi(2,-1);
  for(size_t i=0;i<conns.size();i+=2){
    vi[0] = (conns[i]);
    vi[1] = (conns[i+1]);
    B_PIS->tryInsert(vi);
  }

  // add InteractionGroup to map
  this->m_bpis.insert(make_pair(name,B_PIS));

  console.XDebug()  << "end TSubLattice<T>::addRotBondedIG()\n";
}
 
template <class T>
void TRotSubLattice<T>::addBrittleBeamSCIG()
{
    std::cout << "TRotSubLattice<T>::addBrittleBeamSCIG()\n";

    CVarMPIBuffer param_buffer(this->m_comm);
      
    param_buffer.receiveBroadcast(0);
    int tag=param_buffer.pop_int();
    string name=string(param_buffer.pop_string());
    double kr = param_buffer.pop_double();
    double ks = param_buffer.pop_double();
    double kt = param_buffer.pop_double();
    double kb = param_buffer.pop_double();   
    double cohesion = param_buffer.pop_double();   
    double tCutoff = param_buffer.pop_double();
    double fAngle = param_buffer.pop_double();

    std::cout << "Got prms - kr: " << kr << "  ks: " << ks << " kt: " << kt << " kb: " << kb << " tS: " << cohesion << " tC: " << tCutoff << " fA: " << fAngle << "\n";
    BrittleBeamSCIGP prms(name,kr,ks,kt,kb,cohesion,fAngle,tCutoff,tag);
   
    ParallelInteractionStorage_EB<CRotParticle,BrittleBeamSCInteraction> *B_PIS = 
        new ParallelInteractionStorage_EB<CRotParticle,BrittleBeamSCInteraction>(this->m_ppa,prms);
    
    // construct interactions from connections
    vector<int> conns = TSubLattice<T>::m_temp_conn[tag];
    vector<int> vi(2,-1);
    for(size_t i=0;i<conns.size();i+=2){
        vi[0] = (conns[i]);
        vi[1] = (conns[i+1]);
        B_PIS->tryInsert(vi);
    }
 
    // add InteractionGroup to map
    this->m_bpis.insert(make_pair(name,B_PIS));
    
    std::cout << "end TRotSubLattice<T>::addBrittleBeamSCIG()\n";
}

template <class T>
void TRotSubLattice<T>::addBrittleBeamDZCIG()
{
    std::cout << "TRotSubLattice<T>::addBrittleBeamDZCIG()\n";

    CVarMPIBuffer param_buffer(this->m_comm);

    param_buffer.receiveBroadcast(0);
    int tag=param_buffer.pop_int();
    string name=string(param_buffer.pop_string());
    double kr = param_buffer.pop_double();
    double ks = param_buffer.pop_double();
    double kt = param_buffer.pop_double();
    double kb = param_buffer.pop_double();
    double cohesion = param_buffer.pop_double();
    double tCutoff = param_buffer.pop_double();
    double cCutoff = param_buffer.pop_double();
    double fAngle = param_buffer.pop_double();
    double beta1 = param_buffer.pop_double();
    double beta2 = param_buffer.pop_double();

    std::cout << "Got prms - kr: " << kr << "  ks: " << ks << " kt: " << kt << " kb: " << kb << " tS: " << cohesion << " tC: " << tCutoff << " cC: " << cCutoff << " fA: " << fAngle << " b1: " << beta1 << " b2: " << beta2 << "\n";
    BrittleBeamDZCIGP prms(name,kr,ks,kt,kb,cohesion,fAngle,tCutoff,cCutoff,beta1,beta2,tag);

    ParallelInteractionStorage_EB<CRotParticle,BrittleBeamDZCInteraction> *B_PIS =
        new ParallelInteractionStorage_EB<CRotParticle,BrittleBeamDZCInteraction>(this->m_ppa,prms);

    // construct interactions from connections
    vector<int> conns = TSubLattice<T>::m_temp_conn[tag];
    vector<int> vi(2,-1);
    for(size_t i=0;i<conns.size();i+=2){
        vi[0] = (conns[i]);
        vi[1] = (conns[i+1]);
        B_PIS->tryInsert(vi);
    }

    // add InteractionGroup to map
    this->m_bpis.insert(make_pair(name,B_PIS));

    std::cout << "end TRotSubLattice<T>::addBrittleBeamDZCIG()\n";
}

/*!
  Add thermal bonded interaction group to the lattice. Receive the parameters from master.
  The bonds are created from the neighbor table.
*/
template <class T>
void TRotSubLattice<T>::addRotThermBondedIG()
{
  console.XDebug()  << "TRotSubLattice<T>::addRotThermBondedIG()\n";
  CVarMPIBuffer param_buffer(this->m_comm);
  vector<int> conns;

  // get params
  param_buffer.receiveBroadcast(0);
  int tag=param_buffer.pop_int();
  string name=string(param_buffer.pop_string());
  double kr = param_buffer.pop_double();
  double ks = param_buffer.pop_double();
  double kt = param_buffer.pop_double();
  double kb = param_buffer.pop_double();
  double max_nForce  = param_buffer.pop_double();
  double max_shForce = param_buffer.pop_double();
  double max_tMoment = param_buffer.pop_double();
  double max_bMoment = param_buffer.pop_double();
  double diffusivity = param_buffer.pop_double();

  conns = TSubLattice<T>::m_temp_conn[tag];

  console.XDebug()
    << "Got RotThermBondedIG parameters: tag=" << tag
    << ", name=" << name.c_str()
    << ", kr=" << kr
    << ", ks=" << ks
    << ", kt=" << kt
    << ", kb=" << kb
    << ", nFrc=" << max_nForce
    << ", sFrc=" << max_shForce
    << ", tMom=" << max_tMoment
    << ", bMom=" << max_bMoment
    << ", diffusivity=" << diffusivity
    << "\n";
  // setup InteractionGroup
  CRotThermBondedIGP param;
  param.tag = tag;
  param.kr=kr;
  param.ks=ks;
  param.kt=kt;
  param.kb=kb;
  param.max_nForce = max_nForce ;
  param.max_shForce = max_shForce ;
  param.max_tMoment = max_tMoment ;
  param.max_bMoment = max_bMoment ;
  param.diffusivity = diffusivity ;
  ParallelInteractionStorage_E<CRotThermParticle,CRotThermBondedInteraction> *B_PIS = 
    new ParallelInteractionStorage_EB<CRotThermParticle,CRotThermBondedInteraction>(this->m_ppa,param);

 // recv broadcast connection data
 // console.XDebug() << "rank=" << this->m_tml_comm.rank() << "TSubLattice<T>::addRotThermBondedIG(): receiving conn_data.\n";
 // vector<int> conn_data;
 // this->m_tml_comm.recv_broadcast_cont(conn_data,0);
 // console.XDebug() << "rank=" <<this->m_tml_comm.rank() << "TSubLattice<T>::addRotThermBondedIG(): conn_data.size()=" << conn_data.size() << "\n";
 vector<int> vi(2,-1);
  for(size_t i=0;i<conns.size();i+=2){
    vi[0] = (conns[i]);
    vi[1] = (conns[i+1]);
    B_PIS->tryInsert(vi);
  }

  // add InteractionGroup to map
  this->m_bpis.insert(make_pair(name,B_PIS));

  console.XDebug()  << "end TRotSubLattice<T>::addRotThermBondedIG()\n";
}
