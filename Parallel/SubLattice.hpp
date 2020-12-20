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

// -- project includes -

#include "Parallel/SubLattice.h"
#include "Parallel/MpiInfo.h"
#include "Parallel/mpivbuf.h"
#include "Parallel/mpisgvbuf.h"
#include "Parallel/mpibarrier.h"
#include "Parallel/mpia2abuf.h"
#include "Model/BondedInteraction.h"
#include "Model/CappedBondedInteraction.h"
#include "Model/ShortBondedInteraction.h"
#include "Model/DampingIGP.h"
#include "Model/Damping.h"
#include "Model/RotDamping.h"
#include "Model/ABCDampingIGP.h"
#include "Model/ABCDamping.h"
#include "Model/LocalDampingIGP.h"
#include "Model/LocalDamping.h"
#include "Model/RotLocalDamping.h"
#include "Model/FrictionInteraction.h"
#include "Model/SpringDashpotFrictionInteraction.h"
#include "Model/FractalFriction.h"
#include "Model/AdhesiveFriction.h"
#include "Model/VWFrictionInteraction.h"
#include "Model/RotFricInteraction.h"
#include "Model/RotElasticInteraction.h"
#include "Model/RotThermFricInteraction.h"
#include "Model/RotThermElasticInteraction.h"
#include "Model/RotThermBondedInteraction.h"
#include "Model/HertzianElasticInteraction.h"
#include "Model/HertzianViscoElasticFrictionInteraction.h"
#include "Model/HertzianViscoElasticInteraction.h"
#include "Model/HertzMindlinInteraction.h"
#include "Model/HertzMindlinViscoInteraction.h"
#include "Model/LinearDashpotInteraction.h"
#include "Model/MeshData.h"
#include "Model/MeshData2D.h"
#include "Model/ETriMeshInteraction.h"
#include "Model/BTriMeshInteraction.h"
#include "Model/BMesh2DInteraction.h"
#include "Model/EMesh2DInteraction.h"
#include "Model/TaggedEWallInteractionGroup.h"
#include "Model/ESphereBodyInteractionGroup.h"
#include "Model/FluidCell.h" //fluid contents
#include "Model/FluidInteraction.h" //fluid contents
#include "Parallel/GMRESSolverSlave.h" //fluid contents

// --- parallel storage includes ---
#include "ppa/src/pp_array.h"
#include "pis/pi_storage_eb.h"
#include "pis/pi_storage_ed.h"
#include "pis/pi_storage_ed_t.h"
#include "pis/pi_storage_ne.h"
#include "pis/pi_storage_ne_t.h"
#include "pis/pi_storage_single.h"
#include "pis/pi_storage.h" //fluid contents
#include "pis/trimesh_pis.h"
#include "pis/trimesh_pis_ne.h"
#include "pis/trimesh_pis_eb.h"
#include "pis/mesh2d_pis_eb.h"
#include "pis/mesh2d_pis_ne.h"
#include "Fields/ScalarFluidFieldSlave.h" //fluid contents
#include "Fields/VectorFluidFieldSlave.h" //fluid contents
#include "Fields/ScalarFluidInteractionFieldSlave.h" //fluid contents
#include "Fields/VectorFluidInteractionFieldSlave.h" //fluid contents

// --- field includes ---
#include "Fields/ScalarParticleFieldSlaveTagged.h"
#include "Fields/VectorParticleFieldSlaveTagged.h"
#include "Fields/ScalarInteractionFieldSlaveTagged.h"
#include "Fields/ScalarParticleFieldSlaveTagged.h"
#include "Fields/ScalarInteractionFieldSlaveTagged.h"
#include "Fields/CheckedScalarInteractionFieldSlave.h"
#include "Fields/CheckedScalarInteractionFieldSlaveTagged.h"
#include "Fields/VectorTriangleFieldSlave.h"
#include "Fields/ScalarTriangleFieldSlave.h"
#include "Fields/VectorEdge2DFieldSlave.h"


#include "Model/BodyForceGroup.h"

#include <mpi.h>
#include <stdlib.h>
#include <stdexcept>

// -- STL includes --
#include <algorithm>
#include <stdexcept>
#include <boost/limits.hpp>
using std::runtime_error;

// -- IO includes --
#include <iostream>
using std::cerr;
using std::flush;
using std::endl;
using esys::lsm::CLatticeParam;

//----------------------------------------------
//   TSubLattice functions
//----------------------------------------------

/*!
  construct SubLattice

  \param param Lattice parameters
  \param rank the MPI rank
  \param comm the MPI communicator
*/
template <class T>
TSubLattice<T>::TSubLattice(
  const CLatticeParam &param,
  int rank,
  MPI_Comm comm,
  MPI_Comm worker_comm
)
  : m_ppa(NULL),
    m_dpis(),
    m_bpis(),
    m_singleParticleInteractions(),
    m_damping(),
    m_WIG(),
    m_SIG(),
    m_mesh(),
    m_mesh2d(),
    m_dt(0),
    m_nrange(0),
    m_alpha(0),
    m_last_ns(0),
    m_temp_conn(),
    m_rank(0),
    m_comm(MPI_COMM_NULL),
    m_tml_comm(MPI_COMM_NULL),
    m_worker_comm(MPI_COMM_NULL),
    m_tml_worker_comm(MPI_COMM_NULL),
    m_dims(3, 0),
    packtime(0),
    unpacktime(0),
    commtime(0.0),
    forcetime(0.0),
    m_field_slaves(),
    m_pTimers(NULL)
{
  // cout << "TSubLattice<T>::TSubLattice at " << rank << endl << flush;
  // -- MPI stuff  --
  m_rank=rank;

  // set global communicator
  m_comm=comm;
  m_tml_comm.setComm(m_comm);

  m_dims = param.processDims();

  m_worker_comm=worker_comm;
  //  MPI_Comm_size(m_worker_comm,&m_num_workers);
  m_tml_worker_comm.setComm(m_worker_comm);


  // -- set parameters
  m_alpha=param.alpha();
  m_nrange=param.nrange();
  // cout << "dt,nrange,alpha : " << m_dt << " , " << m_nrange << " , " << m_alpha << "\n";

  commtime=0.0;
  packtime=0.0;
  unpacktime=0.0;
  forcetime=0.0;

  m_last_ns=-1;
}

template <class T>
TSubLattice<T>::~TSubLattice()
{
  console.Debug() << "TSubLattice<T>::~TSubLattice(): enter\n";
  console.Debug()
    << "TSubLattice<T>::~TSubLattice():"
    << " deleting wall interaction groups...\n";
  for(
      typename map<string,AWallInteractionGroup<T>*>::iterator vit=m_WIG.begin();
    vit!=m_WIG.end();
    vit++
  )
  {
    delete vit->second;
  }
  for(
      typename map<string,ASphereBodyInteractionGroup<T>*>::iterator vit=m_SIG.begin();
    vit!=m_SIG.end();
    vit++
  )
  {
    delete vit->second;
  }
  console.Debug()
    << "TSubLattice<T>::~TSubLattice():"
    << " deleting particle array...\n";
  if (m_ppa != NULL) delete m_ppa;
  console.Debug() << "TSubLattice<T>::~TSubLattice(): exit\n";
}

/****fluid contents: begin****/
/*!
  Add fluid interaction to the sublattice. Parameters received from the master
*/
template <class T>
void TSubLattice<T>::addFluidInteraction()
{
  console.XDebug() << "TSubLattice<T>::addFluidInteraction: enter\n" ;
  CVarMPIBuffer param_buffer(m_comm);
  param_buffer.receiveBroadcast(0);

  double cellside=param_buffer.pop_double();
  double Bw=param_buffer.pop_double();
  double Bp=param_buffer.pop_double();
  double Mu=param_buffer.pop_double();
  double alpha=param_buffer.pop_double();
  double flowrate=param_buffer.pop_double();
  double pressure=param_buffer.pop_double();
  Vec3 inflow,outflow;
  inflow[0]=param_buffer.pop_double();
  inflow[1]=param_buffer.pop_double();
  inflow[2]=param_buffer.pop_double();
  outflow[0]=param_buffer.pop_double();
  outflow[1]=param_buffer.pop_double();
  outflow[2]=param_buffer.pop_double();
  m_ppa->setFluid(cellside,Bw,Bp,Mu,alpha,flowrate,pressure,inflow,outflow);

  m_fluidinteraction=new ParallelInteractionStorage_F<T>(m_ppa);
  m_fluidinteraction->update();

  console.XDebug() << "TSubLattice<T>::addFluidInteraction: exit\n" ;
}

/*!
  Add fluid interaction to the sublattice. Parameters received from the master
*/
template <class T>
void TSubLattice<T>::addFluidInteractionVec3()
{
  console.XDebug() << "TSubLattice<T>::addFluidInteractionVec3: enter\n" ;
  CVarMPIBuffer param_buffer(m_comm);
  param_buffer.receiveBroadcast(0);

  Vec3 cellside;
  cellside[0]=param_buffer.pop_double();
  cellside[1]=param_buffer.pop_double();
  cellside[2]=param_buffer.pop_double();
  double Bw=param_buffer.pop_double();
  double Bp=param_buffer.pop_double();
  double Mu=param_buffer.pop_double();
  double alpha=param_buffer.pop_double();
  double flowrate=param_buffer.pop_double();
  double pressure=param_buffer.pop_double();
  Vec3 inflow,outflow;
  inflow[0]=param_buffer.pop_double();
  inflow[1]=param_buffer.pop_double();
  inflow[2]=param_buffer.pop_double();
  outflow[0]=param_buffer.pop_double();
  outflow[1]=param_buffer.pop_double();
  outflow[2]=param_buffer.pop_double();
  m_ppa->setFluidVec3(cellside,Bw,Bp,Mu,alpha,flowrate,pressure,inflow,outflow);

  m_fluidinteraction=new ParallelInteractionStorage_F<T>(m_ppa);
  m_fluidinteraction->update();

  console.XDebug() << "TSubLattice<T>::addFluidInteractionVec3: exit\n" ;
}


/*!
  Update the states of fluid cells
*/

template <class T>
void TSubLattice<T>::updateFluid()
{
  m_ppa->updateFluid();
}

/*!
  Exchange boundary fluid cells
*/

template <class T>
void TSubLattice<T>::exchangeCells()
{
  CVarMPIBuffer buffer(m_comm);
  buffer.receiveBroadcast(0);
  int nt=buffer.pop_int();
  m_ppa->exchangeCells(m_dt,nt);
}


/*!
  add scalar per-fluidcell field to saver list
*/
template <class T>
void TSubLattice<T>::addScalarFluidField()
{
  string fieldname;
  int id;

  m_tml_comm.recv_broadcast_cont(fieldname,0);
  m_tml_comm.recv_broadcast(id,0);

  CFluidCell::ScalarFieldFunction rdf=CFluidCell::getScalarFieldFunction(fieldname);
  ScalarFluidFieldSlave<T> *new_sffs;
  new_sffs=new ScalarFluidFieldSlave<T>(&m_tml_comm,m_ppa,rdf);

  m_field_slaves.insert(make_pair(id,new_sffs));
}

/*!
  add vector per-fluidcell field to saver list
*/
template <class T>
void TSubLattice<T>::addVectorFluidField()
{
  string fieldname;
  int id;

  m_tml_comm.recv_broadcast_cont(fieldname,0);
  m_tml_comm.recv_broadcast(id,0);

  CFluidCell::VectorFieldFunction rdf=CFluidCell::getVectorFieldFunction(fieldname);
  VectorFluidFieldSlave<T> *new_vffs;
  new_vffs=new VectorFluidFieldSlave<T>(&m_tml_comm,m_ppa,rdf);

  m_field_slaves.insert(make_pair(id,new_vffs));
}


/*!
  add scalar fluid interaction field to saver list
*/
template <class T>
void TSubLattice<T>::addScalarFluidInteractionField()
{
  console.XDebug() << "TSubLattice<T>::addScalarFluidInteractionField\n";
  string fieldname;
  int id;

  m_tml_comm.recv_broadcast_cont(fieldname,0);
  console.XDebug() << "recvd. fieldname: " << fieldname << "\n";
  m_tml_comm.recv_broadcast(id,0);
  console.XDebug() << "recvd. id: " << id << "\n";

  CFluidInteraction::ScalarFieldFunction rdf=CFluidInteraction::getScalarFieldFunction(fieldname);
  ScalarFluidInteractionFieldSlave<T> *new_sfifs;
  new_sfifs=new ScalarFluidInteractionFieldSlave<T>(&m_tml_comm,m_fluidinteraction,rdf);
  m_field_slaves.insert(make_pair(id,new_sfifs));

  console.XDebug() << "end TSubLattice<T>::addScalarFluidInteractionField\n";
}

/*!
  add vector fluid interaction field to saver list
*/
template <class T>
void TSubLattice<T>::addVectorFluidInteractionField()
{
  console.XDebug() << "TSubLattice<T>::addVectorFluidInteractionField\n";
  string fieldname;
  int id;

  m_tml_comm.recv_broadcast_cont(fieldname,0);
  console.XDebug() << "recvd. fieldname: " << fieldname << "\n";
  m_tml_comm.recv_broadcast(id,0);
  console.XDebug() << "recvd. id: " << id << "\n";

  CFluidInteraction::VectorFieldFunction rdf=CFluidInteraction::getVectorFieldFunction(fieldname);
  VectorFluidInteractionFieldSlave<T> *new_vfifs;
  new_vfifs=new VectorFluidInteractionFieldSlave<T>(&m_tml_comm,m_fluidinteraction,rdf);
  m_field_slaves.insert(make_pair(id,new_vfifs));

  console.XDebug() << "end TSubLattice<T>::addVectorFluidInteractionField\n";
}

/*!
  send coefficients to master for solving linear equations
*/
template <class T>
void TSubLattice<T>::sendCoeffi()
{
  string fieldname;
  m_tml_comm.recv_broadcast_cont(fieldname,0);

  CFluidCell::ScalarFieldFunction rdf=CFluidCell::getScalarFieldFunction(fieldname);
  vector<pair<Vec3,double> > data_vec;
  data_vec=m_ppa->forAllInnerCellsGetIndexed(rdf);
  m_tml_comm.send_gather(data_vec,0);

  data_vec.clear();
}

/*!
  receive pore pressures from master
*/
template <class T>
void TSubLattice<T>::recvPressure()
{
  vector<pair<Vec3,double> > pressure;
  m_tml_comm.recv_broadcast_cont(pressure,0);
  m_ppa->setPressure(pressure);
  pressure.clear();
}

/*!
  solve matrix by GMRES solver
*/
template <class T>
void TSubLattice<T>::solveMatrix()
{
  vector<int> row,col;
  vector<double> b,v;
  int number,nz;
  m_tml_comm.recv_broadcast_cont(b,0);
  m_tml_comm.recv_broadcast_cont(row,0);
  m_tml_comm.recv_broadcast_cont(col,0);
  m_tml_comm.recv_broadcast_cont(v,0);
  m_tml_comm.recv_broadcast(number,0);
  m_tml_comm.recv_broadcast(nz,0);

  m_solver=new GMRESSolverSlave(&m_worker_comm,&m_tml_comm,b,row,col,v,number,nz);
  m_solver->LocalMatrixSolving();
  //clean up
  row.clear();col.clear();
  b.clear();v.clear();
  delete m_solver;
}

/****fluid contents: end****/


/*!
  Initialize particle storage. The dimensions are global.

  \param min minimum corner
  \param max maximum corner
*/
template <class T>
void TSubLattice<T>::initNeighborTable(const Vec3& min,const Vec3& max)
{
  console.XDebug() << "TSubLattice<T>::initNeighborTable(" << min <<  "," << max << ")\n";
  // make size fit range
  double xsize=max.X()-min.X();
  xsize=m_nrange*ceil(xsize/m_nrange);
  double ysize=max.Y()-min.Y();
  ysize=m_nrange*ceil(ysize/m_nrange);
  double zsize=max.Z()-min.Z();
  zsize=m_nrange*ceil(zsize/m_nrange);
  Vec3 grow=Vec3(xsize,ysize,zsize)-(max-min); // size increase
  Vec3 nmin=min-0.5*grow; // distribute symmetrically
  Vec3 nmax=max+0.5*grow;
  console.XDebug() << "range=" << m_nrange << ", new min,max: " << nmin << ", " << nmax << "\n";

  // construct particle array
  TML_Comm *ntcomm=new TML_Comm(m_worker_comm);
  m_ppa=new  ParallelParticleArray<T>(ntcomm,m_dims,nmin,nmax,m_nrange,m_alpha);
  //m_ppa=new  ParallelParticleArray<T>(ntcomm,3,nmin,nmax,m_nrange);

  //  makeFields(); // put here to make sure ppa is initialized before makeFields

  console.XDebug() << "end TSubLattice<T>::initNeighborTable\n";
}

template <class T>
void TSubLattice<T>::do2dCalculations(bool do2d)
{
  T::setDo2dCalculations(do2d);
}

template <class T>
int TSubLattice<T>::getNumParticles()
{
  return m_ppa->getInnerSize();
}

/*!
  Initialize particle storage with some circular boundaries.
  The dimensions are global.

  \param min minimum corner
  \param max maximum corner
*/
template <class T>
void TSubLattice<T>::initNeighborTable(const Vec3& min,const Vec3& max,const vector<bool>& circ)
{
  console.XDebug() << "TSubLattice<T>::initNeighborTable(" << min <<  "," << max << ") circ\n";
  double xsize,ysize,zsize;
  // if dimension is circular, check if size fits range, otherwise make it fit
  // x - dim
  if(circ[0])
  {
    xsize=max.X()-min.X();
    if(fabs(xsize/m_nrange-lrint(xsize/m_nrange))>1e-6)
    {
      //console.Critical() << "circular x-dimension doesn't fit range !\n";
      console.Info() << "Circular x-size incompatible with range, adjusting...\n";
      m_nrange = xsize/floor(xsize/m_nrange);
      console.Info() << "New range = " << m_nrange << "\n";
    }
    //xsize+=2.0*m_nrange; // padding on the circular ends
  }
  else
  {
    xsize=max.X()-min.X();
    xsize=m_nrange*ceil(xsize/m_nrange);
  }
  // y - dim
  if(circ[1])
  {
    ysize=max.Y()-min.Y();
    if(fabs(ysize/m_nrange-lrint(ysize/m_nrange))>1e-6)
    {
      console.Critical() << "circular y-dimension doesn't fit range !\n";
    }
    ysize+=2.0*m_nrange; // padding on the circular ends
  }
  else
  {
    ysize=max.Y()-min.Y();
    ysize=m_nrange*ceil(ysize/m_nrange);
  }
  // z - dim
  if(circ[2])
  {
    zsize=max.Z()-min.Z();
    if(fabs(zsize/m_nrange-lrint(zsize/m_nrange))>1e-6)
    {
      console.Critical() << "circular z-dimension doesn't fit range !\n";
    }
    zsize+=2.0*m_nrange; // padding on the circular ends
  }
  else
  {
    zsize=max.Z()-min.Z();
    zsize=m_nrange*ceil(zsize/m_nrange);
  }
  Vec3 grow=Vec3(xsize,ysize,zsize)-(max-min); // size increase
  Vec3 nmin=min-0.5*grow; // distribute symmetrically
  Vec3 nmax=max+0.5*grow;
  console.XDebug() << "range, new min, max: " << m_nrange << " " << nmin << nmax << "\n";
  // construct particle array
  TML_Comm *ntcomm=new TML_Comm(m_worker_comm);
  m_ppa=new ParallelParticleArray<T>(ntcomm,m_dims,circ,nmin,nmax,m_nrange,m_alpha);

  //  makeFields(); // put here to make sure ppa is initialized before makeFields

  console.XDebug() << "end TSubLattice<T>::initNeighborTable (circ)\n";
}

/*!
  Receive particle from a TML Communicator

  \param comm the Communicator
  \warning the type of particle is not checked
*/
template <class T>
void TSubLattice<T>::receiveParticles()
{
  //cout<<"begin TSubLattice<T>::receiveParticles() at rank "<<m_rank<<endl;
  console.XDebug() << "TSubLattice<T>::receiveParticles: enter\n";

  vector<T> recv_buffer;
  CMPIBarrier barrier(m_comm);

  m_tml_comm.recv_broadcast_cont_packed(recv_buffer,0);
  console.XDebug() << "recvd " << recv_buffer.size() << " particles \n";
  m_ppa->insert(recv_buffer);

  barrier.wait("TSubLattice<T>::receiveParticles");

  console.XDebug() << "TSubLattice<T>::receiveParticles: exit\n";
  //cout<<"end TSubLattice<T>::receiveParticles() at rank "<<m_rank<<endl;
}


/*!
  Receive connections from a TML Communicator

  \param comm the Communicator
*/
template <class T>
void TSubLattice<T>::receiveConnections()
{
  console.XDebug() << "TSubLattice<T>::receiveConnections: enter\n";

  vector<int> recv_buffer;
  CMPIBarrier barrier(m_comm);

  m_tml_comm.recv_broadcast_cont_packed(recv_buffer,0);
  console.XDebug() << "recvd " << recv_buffer.size() << " connections \n";
  vector<int>::iterator it;
  for (it = recv_buffer.begin(); it != recv_buffer.end(); it+=3)
  {
    if ( (m_ppa->getParticlePtrByIndex( *(it+1)) == NULL ) ||
         (m_ppa->getParticlePtrByIndex( *(it+2)) == NULL ) )
    {
      continue;
    }
    m_temp_conn[*(it)].push_back(*(it+1));
    m_temp_conn[*(it)].push_back(*(it+2));
  }

  barrier.wait("TSubLattice<T>::receiveConnections");

  console.XDebug() << "TSubLattice<T>::receiveConnections: exit\n";
}


/*!
  Add wall to the sublattice. Parameters received from the master
*/
template <class T>
void TSubLattice<T>::addWall()
{
  console.XDebug() << "TSubLattice<T>::addWall: enter\n" ;
  CVarMPIBuffer param_buffer(m_comm);
  param_buffer.receiveBroadcast(0);

  string name=param_buffer.pop_string();
  Vec3 ipos=param_buffer.pop_vector();
  Vec3 inorm=param_buffer.pop_vector();

  m_walls[name]=new CWall(ipos,inorm);

  console.XDebug() << "TSubLattice<T>::addWall: exit\n" ;
}

/*!
  Add sphere body to the sublattice. Parameters received from the master
*/
template <class T>
void TSubLattice<T>::addSphereBody()
{
  console.XDebug() << "TSubLattice<T>::addSphereBody: enter\n" ;
  CVarMPIBuffer param_buffer(m_comm);
  param_buffer.receiveBroadcast(0);

  string name=param_buffer.pop_string();
  Vec3 ipos=param_buffer.pop_vector();
  double radius=param_buffer.pop_double();

  m_spheres[name]=new CSphereBody(ipos,radius);

  console.XDebug() << "TSubLattice<T>::addSphereBody: exit\n" ;
}

/*!
  add elastic wall interaction group. Parameters received from the master
*/
template <class T>
void TSubLattice<T>::addElasticWIG()
{
  console.XDebug() << "TSubLattice<T>::addElasticWIG: enter\n" ;
  CVarMPIBuffer param_buffer(m_comm);

  // receive params from master
  param_buffer.receiveBroadcast(0);

  CEWallIGP* wigp=extractEWallIGP(&param_buffer);

  string wallname=wigp->getWallName();
  map<string,CWall*>::iterator iter=m_walls.find(wallname);
  if(iter!=m_walls.end()){
    AWallInteractionGroup<T>* newCEWIG =
      new CEWallInteractionGroup<T>(
        &m_tml_worker_comm,
        m_walls[wallname],
        wigp
      );
    newCEWIG->Update(m_ppa);
    m_WIG.insert(make_pair(wigp->getName(),newCEWIG));
  } else {
    std::stringstream msg;
    msg << "wall name '" << wallname << "' not found in map of walls";
    throw std::runtime_error(msg.str().c_str());
  }

  delete wigp;
  console.XDebug() << "TSubLattice<T>::addElasticWIG: exit\n" ;
}

/*!
  add elastic sphere body interaction group. Parameters received from the master
*/
template <class T>
void TSubLattice<T>::addESphereBodyIG()
{
  console.XDebug() << "TSubLattice<T>::addESphereBodyIG: enter\n" ;
  CVarMPIBuffer param_buffer(m_comm);

  // receive params from master
  param_buffer.receiveBroadcast(0);

  CESphereBodyIGP* wigp=extractESphereBodyIGP(&param_buffer);

  string wallname=wigp->getSphereBodyName();
  map<string,CSphereBody*>::iterator iter=m_spheres.find(wallname);
  if(iter!=m_spheres.end()){
    ASphereBodyInteractionGroup<T>* newCEWIG =
      new CESphereBodyInteractionGroup<T>(
        &m_tml_worker_comm,
        m_spheres[wallname],
        wigp
      );
    newCEWIG->Update(m_ppa);
    m_SIG.insert(make_pair(wigp->getName(),newCEWIG));
  } else {
    std::stringstream msg;
    msg << "sphere body name '" << wallname << "' not found in map of sphere bodies";
    throw std::runtime_error(msg.str().c_str());
  }

  delete wigp;
  console.XDebug() << "TSubLattice<T>::addESphereBodyIG: exit\n" ;
}

/*!
  add tagged elastic wall interaction group. Parameters received from the master
*/
template <class T>
void TSubLattice<T>::addTaggedElasticWIG()
{
  console.XDebug() << "TSubLattice<T>::addTaggedElasticWIG: enter\n" ;
  CVarMPIBuffer param_buffer(m_comm);

  // receive params from master
  param_buffer.receiveBroadcast(0);

  CEWallIGP* wigp=extractEWallIGP(&param_buffer);
  int tag=param_buffer.pop_int();
  int mask=param_buffer.pop_int();

  string wallname=wigp->getWallName();

  console.XDebug() << wallname << " tag= " << tag << " mask= " << mask <<"\n"	;
  map<string,CWall*>::iterator iter=m_walls.find(wallname);
  if(iter!=m_walls.end()){
     AWallInteractionGroup<T>* newCTEWIG =
	new CTaggedEWallInteractionGroup<T>(
		&m_tml_worker_comm,
		m_walls[wallname],
		wigp,
		tag,
		mask
	);
     newCTEWIG->Update(m_ppa);
     m_WIG.insert(make_pair(wigp->getName(),newCTEWIG));
   } else {
     std::stringstream msg;
     msg << "wall name '" << wallname << "' not found in map of walls";
     throw std::runtime_error(msg.str().c_str());
   }

  delete wigp;
  console.XDebug() << "TSubLattice<T>::addTaggedElasticWIG: exit\n" ;
}


/*!
  add bonded wall interaction group. Parameters received from the master
*/
template <class T>
void TSubLattice<T>::addBondedWIG()
{
  console.XDebug() << "TSubLattice<T>::addBondedWIG: enter\n" ;
  CVarMPIBuffer param_buffer(m_comm);

  // receive params from master
  param_buffer.receiveBroadcast(0);

  CBWallIGP* wigp=extractBWallIGP(&param_buffer);

  string wallname=wigp->getWallName();
  map<string,CWall*>::iterator iter=m_walls.find(wallname);
  if(iter!=m_walls.end()){
    AWallInteractionGroup<T>* newCBWIG=new CBWallInteractionGroup<T>(&m_tml_worker_comm,m_walls[wallname],wigp);
    newCBWIG->Update(m_ppa);
    m_WIG.insert(make_pair(wigp->getName(),newCBWIG));
  } else {
    console.Error() << "wall name " << wallname << " not found in map of walls\n";
  }

  delete wigp;
  console.XDebug() << "TSubLattice<T>::addBondedWIG: exit\n" ;
}

/*!
  add bonded wall interaction group with direction dependend elasticity . Parameters received from the master
*/
template <class T>
void TSubLattice<T>::addDirBondedWIG()
{
  console.XDebug() << "TSubLattice<T>::addDirBondedWIG: enter\n" ;
  CVarMPIBuffer param_buffer(m_comm);

  // receive params from master
  param_buffer.receiveBroadcast(0);

  CSoftBWallIGP* wigp=extractSoftBWallIGP(&param_buffer);

  string wallname=wigp->getWallName();
  map<string,CWall*>::iterator iter=m_walls.find(wallname);
  if(iter!=m_walls.end()){
    AWallInteractionGroup<T>* newCDWIG=new CSoftBWallInteractionGroup<T>(&m_tml_worker_comm,m_walls[wallname],wigp);
    newCDWIG->Update(m_ppa);
    m_WIG.insert(make_pair(wigp->getName(),newCDWIG));
  } else {
    console.Error() << "wall name " << wallname << " not found in map of walls\n";
  }

  delete wigp;
  console.XDebug() << "TSubLattice<T>::addDirBondedWIG: exit\n" ;
}

/*!
  get position of a wall
*/
template <class T>
void TSubLattice<T>::getWallPos()
{
  console.XDebug() << "TSubLattice<T>::getWallPosition: enter\n" ;
  CVarMPIBuffer param_buffer(m_comm);
  Vec3 pos;

  // receive params from master
  param_buffer.receiveBroadcast(0);

  std::string wname=param_buffer.pop_string();
  console.XDebug() << "Wall name: " << wname << "\n";

  // find wall
  map<string,CWall*>::iterator iter=m_walls.find(wname);
  if(iter!=m_walls.end()){
    pos=(iter->second)->getPos();
    console.XDebug() << "Wall pos: " << pos << "\n";
  } else {
    pos=Vec3(0.0,0.0,0.0);
  }

  vector<Vec3> vpos;
  vpos.push_back(pos);
  m_tml_comm.send_gather(vpos,0);
  console.XDebug() << "TSubLattice<T>::getWallPosition: exit\n" ;
}

/*!
  get position of a sphere body
*/
template <class T>
void TSubLattice<T>::getSphereBodyPos()
{
  console.XDebug() << "TSubLattice<T>::getSphereBodyPosition: enter\n" ;
  CVarMPIBuffer param_buffer(m_comm);
  Vec3 pos;

  // receive params from master
  param_buffer.receiveBroadcast(0);

  std::string wname=param_buffer.pop_string();
  console.XDebug() << "Sphere name: " << wname << "\n";

  // find sphere
  map<string,CSphereBody*>::iterator iter=m_spheres.find(wname);
  if(iter!=m_spheres.end()){
    pos=(iter->second)->getPos();
    console.XDebug() << "Sphere pos: " << pos << "\n";
  } else {
    pos=Vec3(0.0,0.0,0.0);
  }

  vector<Vec3> vpos;
  vpos.push_back(pos);
  m_tml_comm.send_gather(vpos,0);
  console.XDebug() << "TSubLattice<T>::getSphereBodyPosition: exit\n" ;
}

/*!
  get force acting on a wall
*/
template <class T>
void TSubLattice<T>::getWallForce()
{
  console.XDebug() << "TSubLattice<T>::getWallForce: enter\n" ;
  CVarMPIBuffer param_buffer(m_comm);
  Vec3 force;

  // receive params from master
  param_buffer.receiveBroadcast(0);

  std::string wname=param_buffer.pop_string();
  console.XDebug() << "Wall name: " << wname << "\n";

  // find wall
  map<string,CWall*>::iterator iter=m_walls.find(wname);
  if(iter!=m_walls.end()){
    force=(iter->second)->getForce();
    console.XDebug() << "Wall force: " << force << "\n";
  } else {
    force=Vec3(0.0,0.0,0.0);
  }

  vector<Vec3> vforce;
  vforce.push_back(force);
  m_tml_comm.send_gather(vforce,0);
  console.XDebug() << "TSubLattice<T>::getWallForce: exit\n" ;
}

/*!
  get force acting on a sphere body
*/
template <class T>
void TSubLattice<T>::getSphereBodyForce()
{
  console.XDebug() << "TSubLattice<T>::getSphereBodyForce: enter\n" ;
  CVarMPIBuffer param_buffer(m_comm);
  Vec3 force;

  // receive params from master
  param_buffer.receiveBroadcast(0);

  std::string wname=param_buffer.pop_string();
  console.XDebug() << "Sphere name: " << wname << "\n";

  // find Sphere
  map<string,CSphereBody*>::iterator iter=m_spheres.find(wname);
  if(iter!=m_spheres.end()){
    force=(iter->second)->getForce();
    console.XDebug() << "Sphere force: " << force << "\n";
  } else {
    force=Vec3(0.0,0.0,0.0);
  }

  vector<Vec3> vforce;
  vforce.push_back(force);
  m_tml_comm.send_gather(vforce,0);
  console.XDebug() << "TSubLattice<T>::getSphereBodyForce: exit\n" ;
}

/*!
  add wall with viscous drag
*/
template <class T>
void TSubLattice<T>::addViscWIG()
{
  console.XDebug() << "TSubLattice<T>::addViscWIG: enter\n" ;
  CVarMPIBuffer param_buffer(m_comm);

  // receive params from master
  param_buffer.receiveBroadcast(0);

  CVWallIGP* wigp=extractVWallIGP(&param_buffer);

  string wallname=wigp->getWallName();
  map<string,CWall*>::iterator iter=m_walls.find(wallname);
  if(iter!=m_walls.end()){
    AWallInteractionGroup<T>* newCVWIG=new CViscWallIG<T>(&m_tml_worker_comm,m_walls[wallname],wigp);
    newCVWIG->Update(m_ppa);
    m_WIG.insert(make_pair(wigp->getName(),newCVWIG));
  } else {
    console.Error() << "wall name " << wallname << " not found in map of walls\n";
  }

  delete wigp;
  console.XDebug() << "TSubLattice<T>::addViscWIG: exit\n" ;
}

/*!
  Add a PairInteractionGroup to the lattice
*/
template <class T>
void TSubLattice<T>::addPairIG()
{
  console.XDebug()  << "TSubLattice<T>::addPairIG()\n";
  CVarMPIBuffer param_buffer(m_comm,2000);

  // get params
  param_buffer.receiveBroadcast(0);
  string type = param_buffer.pop_string();
  console.XDebug()<< "PIG type: " << type.c_str() << "\n";
  string name = param_buffer.pop_string();
  console.XDebug()<< "PIG name: " << name.c_str() << "\n";

  doAddPIG(name,type,param_buffer,false);

  console.XDebug()  << "end TSubLattice<T>::addPairIG()\n";
}

/*!
  Add a tagged PairInteractionGroup to the lattice
*/
template <class T>
void TSubLattice<T>::addTaggedPairIG()
{
  console.XDebug()  << "TSubLattice<T>::addTaggedPairIG()\n";
  CVarMPIBuffer param_buffer(m_comm,2000);

  // get params
  param_buffer.receiveBroadcast(0);
  string type = param_buffer.pop_string();
  console.XDebug()<< "PIG type: " << type.c_str() << "\n";
  string name = param_buffer.pop_string();
  console.XDebug()<< "PIG name: " << name.c_str() << "\n";

  doAddPIG(name,type,param_buffer,true);

  console.XDebug()  << "end TSubLattice<T>::addTaggedPairIG()\n";
}

/*!
  do the actual work adding the PIG

  \param name the name of the PIG
  \param type the type of the PIG
  \param param_buffer the buffer containing the rest of the parameters
*/
template <class T>
bool TSubLattice<T>::doAddPIG(const string& name,const string& type,CVarMPIBuffer& param_buffer, bool tagged)
{
  bool res=false;
  AParallelInteractionStorage* new_pis = NULL;

  if(type=="Elastic") {
    CElasticIGP eigp;
    eigp.m_k=param_buffer.pop_double();
    eigp.m_scaling=static_cast<bool>(param_buffer.pop_int());
    if(tagged){
      int tag1=param_buffer.pop_int();
      int mask1=param_buffer.pop_int();
      int tag2=param_buffer.pop_int();
      int mask2=param_buffer.pop_int();
      console.XDebug() << "tag1, mask1, tag2, mask2 "
		       << tag1 << " , " << mask1 << " , "
		       << tag2 << " , " << mask2 << "\n";
      new_pis=new ParallelInteractionStorage_NE_T<T,CElasticInteraction>(m_ppa,eigp,tag1,mask1,tag2,mask2);
    }else{
      new_pis=new ParallelInteractionStorage_NE<T,CElasticInteraction>(m_ppa,eigp);
    }
    res=true;
  } else if (type=="Friction") {
    CFrictionIGP figp;
    figp.k=param_buffer.pop_double();
    figp.mu=param_buffer.pop_double();
    figp.k_s=param_buffer.pop_double();
    figp.dt=param_buffer.pop_double();
    figp.m_scaling=static_cast<bool>(param_buffer.pop_int());
    console.XDebug() << "k,mu,k_s,dt: " << figp.k << " , " << figp.mu << " , "
    << figp.k_s << " , " << figp.dt << "\n";
    if(tagged){
	int tag1=param_buffer.pop_int();
	int mask1=param_buffer.pop_int();
	int tag2=param_buffer.pop_int();
	int mask2=param_buffer.pop_int();
	console.XDebug() << "tag1, mask1, tag2, mask2 " 
			 << tag1 << " , " << mask1 << " , " 
			 << tag2 << " , " << mask2 << "\n";  
	new_pis=new ParallelInteractionStorage_ED_T<T,CFrictionInteraction>(m_ppa,figp,tag1,mask1,tag2,mask2);
    } else {
	new_pis=new ParallelInteractionStorage_ED<T,CFrictionInteraction>(m_ppa,figp);
    }
    res=true;
  } else if (type=="Spring_Dashpot_Friction") {
    CSpringDashpotFrictionIGP figp;
    figp.youngModulus=param_buffer.pop_double();
    figp.poissonRatio=param_buffer.pop_double();
    figp.cor=param_buffer.pop_double();
    figp.mu=param_buffer.pop_double();
    figp.dt=param_buffer.pop_double();
    console.XDebug() << "Y,nu,cor,mu,dt: " << figp.youngModulus << " , " << figp.poissonRatio << " , " << figp.cor << " , "
    << figp.mu << " , " << figp.dt << "\n";
    if(tagged){
        int tag1=param_buffer.pop_int();
        int mask1=param_buffer.pop_int();
        int tag2=param_buffer.pop_int();
        int mask2=param_buffer.pop_int();
        console.XDebug() << "tag1, mask1, tag2, mask2 " 
                         << tag1 << " , " << mask1 << " , "
                         << tag2 << " , " << mask2 << "\n";  
        new_pis=new ParallelInteractionStorage_ED_T<T,CSpringDashpotFrictionInteraction>(m_ppa,figp,tag1,mask1,tag2,mask2);
    } else {
        new_pis=new ParallelInteractionStorage_ED<T,CSpringDashpotFrictionInteraction>(m_ppa,figp);
    }
    res=true;
  } else if (type=="AdhesiveFriction") {
    CAdhesiveFrictionIGP figp;
    figp.k=param_buffer.pop_double();
    figp.mu=param_buffer.pop_double();
    figp.k_s=param_buffer.pop_double();
    figp.dt=param_buffer.pop_double();
    figp.r_cut=param_buffer.pop_double();
    console.XDebug()
      << "k,mu,k_s,dt,r_cut: " << figp.k << " , " << figp.mu << " , "
      << figp.k_s << " , " << figp.dt << " " << figp.r_cut << "\n";
    if(tagged){
	int tag1=param_buffer.pop_int();
	int mask1=param_buffer.pop_int();
	int tag2=param_buffer.pop_int();
	int mask2=param_buffer.pop_int();
	console.XDebug() << "tag1, mask1, tag2, mask2 " 
			 << tag1 << " , " << mask1 << " , " 
			 << tag2 << " , " << mask2 << "\n";  
	new_pis=new ParallelInteractionStorage_ED_T<T,CAdhesiveFriction>(m_ppa,figp,tag1,mask1,tag2,mask2);    
    } else {
	new_pis=new ParallelInteractionStorage_ED<T,CAdhesiveFriction>(m_ppa,figp);
    }
    res=true;
  } else if (type=="FractalFriction") {
    FractalFrictionIGP figp;
    figp.k=param_buffer.pop_double();
    figp.mu_0=param_buffer.pop_double();
    figp.k_s=param_buffer.pop_double();
    figp.dt=param_buffer.pop_double();
    console.XDebug() << "k,mu_0,k_s,dt: " << figp.k << " , " << figp.mu_0 << " , "
                     << figp.k_s << " , " << figp.dt << "\n";
    figp.x0=param_buffer.pop_double();
    figp.y0=param_buffer.pop_double();
    figp.dx=param_buffer.pop_double();
    figp.dy=param_buffer.pop_double();
    figp.nx=param_buffer.pop_int();
    figp.ny=param_buffer.pop_int();
    console.XDebug()
      <<"x0,y0,dx,dy,nx,ny: "
      << figp.x0 << " , " << figp.y0 << " , "
      << figp.dx << " , " << figp.dy << " ,"
      << figp.nx << " , " << figp.ny << "\n";
    figp.mu = boost::shared_ptr<double>(new double[figp.nx*figp.ny]);

    for(int i=0;i<figp.nx*figp.ny;i++)
    {
      (figp.mu.get())[i]=param_buffer.pop_double();
      // console.XDebug() << i << " , " << figp.mu[i] << "\n";
    }
    if(tagged){
	int tag1=param_buffer.pop_int();
	int mask1=param_buffer.pop_int();
	int tag2=param_buffer.pop_int();
	int mask2=param_buffer.pop_int();
	console.XDebug() << "tag1, mask1, tag2, mask2 " 
			 << tag1 << " , " << mask1 << " , " 
			 << tag2 << " , " << mask2 << "\n";  
	new_pis = new ParallelInteractionStorage_ED_T<T,CFractalFriction>(m_ppa,figp,tag1,mask1,tag2,mask2);
    } else {
	new_pis = new ParallelInteractionStorage_ED<T,CFractalFriction>(m_ppa,figp);
    }
    res=true;
 } else if(type=="VWFriction") {
    VWFrictionIGP figp;

    figp.k=param_buffer.pop_double();
    figp.mu=param_buffer.pop_double();
    figp.k_s=param_buffer.pop_double();
    figp.dt=param_buffer.pop_double();
    figp.m_alpha=param_buffer.pop_double();
    console.XDebug()
      << "k,mu,k_s,dt,alpha: " << figp.k << " , " << figp.mu << " , "
      << figp.k_s << " , " << figp.dt << "\n";
    if(tagged){
	int tag1=param_buffer.pop_int();
	int mask1=param_buffer.pop_int();
	int tag2=param_buffer.pop_int();
	int mask2=param_buffer.pop_int();
	console.XDebug() << "tag1, mask1, tag2, mask2 " 
			 << tag1 << " , " << mask1 << " , " 
			 << tag2 << " , " << mask2 << "\n";  
	new_pis=new ParallelInteractionStorage_ED_T<T,CVWFriction>(m_ppa,figp,tag1,mask1,tag2,mask2);
    } else {
	new_pis=new ParallelInteractionStorage_ED<T,CVWFriction>(m_ppa,figp);
    }
    res=true;
  }  else if(type=="RotElastic"){
    CRotElasticIGP reigp;
    reigp.m_kr=param_buffer.pop_double();
    if(tagged){
	int tag1=param_buffer.pop_int();
	int mask1=param_buffer.pop_int();
	int tag2=param_buffer.pop_int();
	int mask2=param_buffer.pop_int();
	console.XDebug() << "tag1, mask1, tag2, mask2 " 
			 << tag1 << " , " << mask1 << " , " 
			 << tag2 << " , " << mask2 << "\n";  
    
	new_pis=new ParallelInteractionStorage_NE_T<CRotParticle,CRotElasticInteraction>(m_ppa,reigp,tag1,mask1,tag2,mask2);
    } else {
	new_pis=new ParallelInteractionStorage_NE<CRotParticle,CRotElasticInteraction>(m_ppa,reigp);
    }
    res=true;
  } else if (type=="RotFriction"){
    CRotFrictionIGP rfigp;
    rfigp.k=param_buffer.pop_double();
    rfigp.mu_s=param_buffer.pop_double();
    rfigp.mu_d=param_buffer.pop_double();
    rfigp.k_s=param_buffer.pop_double();
    rfigp.dt=param_buffer.pop_double();
    rfigp.scaling=static_cast<bool>(param_buffer.pop_int());
    rfigp.meanR_scaling=static_cast<bool>(param_buffer.pop_int());
    console.XDebug()
      << "k,mu_s,mu_d,k_s,dt,scaling: " << rfigp.k << " , "
      << rfigp.mu_s << " , " << rfigp.mu_d << " , "
      << rfigp.k_s << " , " << rfigp.dt << " , " << rfigp.scaling << "\n";
    if(tagged){
      int tag1=param_buffer.pop_int();
      int mask1=param_buffer.pop_int();
      int tag2=param_buffer.pop_int();
      int mask2=param_buffer.pop_int();
      console.XDebug() << "tag1, mask1, tag2, mask2 "
        << tag1 << " , " << mask1 << " , "
        << tag2 << " , " << mask2 << "\n";
      new_pis=new ParallelInteractionStorage_ED_T<CRotParticle,CRotFrictionInteraction>(m_ppa,rfigp,tag1,mask1,tag2,mask2);
    } else {
      new_pis=new ParallelInteractionStorage_ED<CRotParticle,CRotFrictionInteraction>(m_ppa,rfigp);
    }
    res=true;
  } else if (type == "RotThermElastic") {
    CRotThermElasticIGP eigp;
    eigp.m_kr        = param_buffer.pop_double();
    eigp.diffusivity = param_buffer.pop_double();
    console.XDebug()
      << "k=" << eigp.m_kr << " , "
      << "diffusivity=" << eigp.diffusivity << "\n";
    if(tagged){
	int tag1=param_buffer.pop_int();
	int mask1=param_buffer.pop_int();
	int tag2=param_buffer.pop_int();
	int mask2=param_buffer.pop_int();
	console.XDebug() << "tag1, mask1, tag2, mask2 " 
			 << tag1 << " , " << mask1 << " , " 
			 << tag2 << " , " << mask2 << "\n";  
	    
	new_pis=new ParallelInteractionStorage_NE_T<CRotThermParticle,CRotThermElasticInteraction>(m_ppa,eigp,tag1,mask1,tag2,mask2);
    } else {
	new_pis=new ParallelInteractionStorage_NE<CRotThermParticle,CRotThermElasticInteraction>(m_ppa,eigp);
    }    
    res=true;
  } else if (type == "RotThermFriction") {
    CRotThermFrictionIGP rfigp;
    rfigp.k=param_buffer.pop_double();
    rfigp.mu_s=param_buffer.pop_double();
    rfigp.mu_d=param_buffer.pop_double();
    rfigp.k_s=param_buffer.pop_double();
    rfigp.diffusivity=param_buffer.pop_double();
    rfigp.dt=param_buffer.pop_double();
    console.XDebug()
      << "k=" << rfigp.k << " , "
      << "mu_d=" << rfigp.mu_d << " , "
      << "mu_s=" << rfigp.mu_s << " , "
      << "k_s=" << rfigp.k_s << " , "
      << "diffusivity=" << rfigp.diffusivity << " , "
      << "dt=" << rfigp.dt << "\n";
    if(tagged){
	int tag1=param_buffer.pop_int();
	int mask1=param_buffer.pop_int();
	int tag2=param_buffer.pop_int();
	int mask2=param_buffer.pop_int();
	console.XDebug() << "tag1, mask1, tag2, mask2 " 
			 << tag1 << " , " << mask1 << " , " 
			 << tag2 << " , " << mask2 << "\n";  
	new_pis=new ParallelInteractionStorage_ED_T<CRotThermParticle,CRotThermFrictionInteraction>(m_ppa,rfigp,tag1,mask1,tag2,mask2);    
    } else {
	new_pis=new ParallelInteractionStorage_ED<CRotThermParticle,CRotThermFrictionInteraction>(m_ppa,rfigp);
    }
    res=true;
 } else if(type=="HertzianElastic") {
    CHertzianElasticIGP heigp;
    heigp.m_E=param_buffer.pop_double();
    heigp.m_nu=param_buffer.pop_double();
    if(tagged){
      int tag1=param_buffer.pop_int();
      int mask1=param_buffer.pop_int();
      int tag2=param_buffer.pop_int();
      int mask2=param_buffer.pop_int();
      new_pis=new ParallelInteractionStorage_NE_T<T,CHertzianElasticInteraction>(m_ppa,heigp,tag1,mask1,tag2,mask2);
    }else{
      new_pis=new ParallelInteractionStorage_NE<T,CHertzianElasticInteraction>(m_ppa,heigp);
    }
    res=true;
  } else if(type=="HertzianViscoElasticFriction") {
    CHertzianViscoElasticFrictionIGP hvefigp;
    hvefigp.m_A=param_buffer.pop_double();
    hvefigp.m_E=param_buffer.pop_double();
    hvefigp.m_nu=param_buffer.pop_double();
    hvefigp.mu=param_buffer.pop_double();
    hvefigp.k_s=param_buffer.pop_double();
    hvefigp.dt=param_buffer.pop_double();
    if(tagged){
      int tag1=param_buffer.pop_int();
      int mask1=param_buffer.pop_int();
      int tag2=param_buffer.pop_int();
      int mask2=param_buffer.pop_int();
      new_pis=new ParallelInteractionStorage_ED_T<T,CHertzianViscoElasticFrictionInteraction>(m_ppa,hvefigp,tag1,mask1,tag2,mask2);
    }else{
      new_pis=new ParallelInteractionStorage_ED<T,CHertzianViscoElasticFrictionInteraction>(m_ppa,hvefigp);
    }
    res=true;
  } else if(type=="HertzianViscoElastic") {
    CHertzianViscoElasticIGP hveigp;
    hveigp.m_A=param_buffer.pop_double();
    hveigp.m_E=param_buffer.pop_double();
    hveigp.m_nu=param_buffer.pop_double();
    if(tagged){
      int tag1=param_buffer.pop_int();
      int mask1=param_buffer.pop_int();
      int tag2=param_buffer.pop_int();
      int mask2=param_buffer.pop_int();
      new_pis=new ParallelInteractionStorage_NE_T<T,CHertzianViscoElasticInteraction>(m_ppa,hveigp,tag1,mask1,tag2,mask2);
    }else{
      new_pis=new ParallelInteractionStorage_NE<T,CHertzianViscoElasticInteraction>(m_ppa,hveigp);
    }
    res=true;
  } else if(type=="HertzMindlin") {
    CHertzMindlinIGP hmigp;
    hmigp.m_E=param_buffer.pop_double();
    hmigp.m_nu=param_buffer.pop_double();
    hmigp.mu=param_buffer.pop_double();
    hmigp.dt=param_buffer.pop_double();
    if(tagged){
      int tag1=param_buffer.pop_int();
      int mask1=param_buffer.pop_int();
      int tag2=param_buffer.pop_int();
      int mask2=param_buffer.pop_int();
      new_pis=new ParallelInteractionStorage_NE_T<T,CHertzMindlinInteraction>(m_ppa,hmigp,tag1,mask1,tag2,mask2);
    }else{
      new_pis=new ParallelInteractionStorage_NE<T,CHertzMindlinInteraction>(m_ppa,hmigp);
    }
    res=true;
  } else if(type=="HertzMindlinVisco") {
    CHertzMindlinViscoIGP hmvigp;
    hmvigp.m_E=param_buffer.pop_double();
    hmvigp.m_nu=param_buffer.pop_double();
    hmvigp.mu=param_buffer.pop_double();
    hmvigp.m_COR=param_buffer.pop_double();
    hmvigp.dt=param_buffer.pop_double();
    if(tagged){
      int tag1=param_buffer.pop_int();
      int mask1=param_buffer.pop_int();
      int tag2=param_buffer.pop_int();
      int mask2=param_buffer.pop_int();
      new_pis=new ParallelInteractionStorage_NE_T<T,CHertzMindlinViscoInteraction>(m_ppa,hmvigp,tag1,mask1,tag2,mask2);
    }else{
      new_pis=new ParallelInteractionStorage_NE<T,CHertzMindlinViscoInteraction>(m_ppa,hmvigp);
    }
    res=true;
  } else if(type=="LinearDashpot") {
    CLinearDashpotIGP heigp;
    heigp.m_damp=param_buffer.pop_double();
    heigp.m_cutoff=param_buffer.pop_double();
    if(tagged){
      int tag1=param_buffer.pop_int();
      int mask1=param_buffer.pop_int();
      int tag2=param_buffer.pop_int();
      int mask2=param_buffer.pop_int();
      new_pis=new ParallelInteractionStorage_NE_T<T,CLinearDashpotInteraction>(m_ppa,heigp,tag1,mask1,tag2,mask2);
    }else{
      new_pis=new ParallelInteractionStorage_NE<T,CLinearDashpotInteraction>(m_ppa,heigp);
    }
    res=true;
  } else {
    cerr  << "Unknown interaction group name "
	  << type
	  << " in TSubLattice<T>::addPairIG()" << endl;
  }

  // add InteractionGroup to map
  if(res) m_dpis.insert(make_pair(name,new_pis));

  return res;
}


/*!
  Add a Triangle Mesh
*/
template <class T>
void TSubLattice<T>::addTriMesh()
{
  console.XDebug()  << "TSubLattice<T>::addTriMesh()\n";

  // MPI buffer
  CVarMPIBuffer param_buffer(m_comm);
  // data buffers
  vector<MeshNodeData> node_recv_buffer;
  vector<MeshTriData> tri_recv_buffer;

  // receive params
  param_buffer.receiveBroadcast(0);

  // extract name
  string name = param_buffer.pop_string();
  console.XDebug()<< "TriMesh name: " << name.c_str() << "\n";

  // receive nodes
  m_tml_comm.recv_broadcast_cont_packed(node_recv_buffer,0);
  console.XDebug() << "recvd " << node_recv_buffer.size() << " nodes \n";

  // receive triangles
  m_tml_comm.recv_broadcast_cont_packed(tri_recv_buffer,0);
  console.XDebug() << "recvd " << tri_recv_buffer.size() << " triangles \n";

  // load mesh into new TriMesh
  TriMesh* new_tm=new TriMesh();
  new_tm->LoadMesh(node_recv_buffer,tri_recv_buffer);

  m_mesh.insert(make_pair(name,new_tm));

  console.XDebug()  << "end TSubLattice<T>::addTriMesh()\n";
}

/*!
  Add a TriMesh interaction group
*/
template <class T>
void TSubLattice<T>::addTriMeshIG()
{
  console.XDebug()  << "TSubLattice<T>::addTriMeshIG()\n";

  // MPI buffer
  CVarMPIBuffer param_buffer(m_comm);

  // receive params
  param_buffer.receiveBroadcast(0);

  // extract name & type
  string type = param_buffer.pop_string();
  console.XDebug()<< "TriMeshIG type: " << type.c_str() << "\n";
  string name = param_buffer.pop_string();
  console.XDebug()<< "TriMeshIG name: " << name.c_str() << "\n";
  string meshname = param_buffer.pop_string();
  console.XDebug()<< "TriMeshIG mesh name: " << meshname.c_str() << "\n";

  // get pointer to mesh
  TriMesh* tmp=NULL;
  if (m_mesh.find(meshname) != m_mesh.end())
  {
    tmp = m_mesh[meshname];
  }
  if(tmp==NULL){
    throw runtime_error("unknown mesh name in TSubLattice<T>::addTriMeshIG:" + meshname);
  }
  // switch on type,extract params & construc new TriMeshIG
  if(type=="Elastic")
  {
    AParallelInteractionStorage* new_pis;
    ETriMeshIP tmi;
    tmi.k=param_buffer.pop_double();
    new_pis = new TriMesh_PIS_NE<T,ETriMeshInteraction>(tmp,m_ppa,tmi);
    m_dpis.insert(make_pair(name,new_pis));
  } else { // unknown type-> throw
    throw runtime_error("unknown type in TSubLattice<T>::addTriMeshIG:" + type);
  }


  console.XDebug()  << "end TSubLattice<T>::addTriMeshIG()\n";
}

/*!
  Add a bonded TriMesh interaction group
*/
template <class T>
void TSubLattice<T>::addBondedTriMeshIG()
{
  console.XDebug()  << "TSubLattice<T>::addBondedTriMeshIG()\n";

  // MPI buffer
  CVarMPIBuffer param_buffer(m_comm);

  // receive params
  BTriMeshIP param;
  param_buffer.receiveBroadcast(0);

  // extract name.meshname & params
  string name = param_buffer.pop_string();
  console.XDebug()<< "BTriMeshIG name: " << name.c_str() << "\n";
  string meshname = param_buffer.pop_string();
  console.XDebug()<< "BTriMeshIG mesh name: " << meshname.c_str() << "\n";
  param.k = param_buffer.pop_double();
  console.XDebug()<< "BTriMeshIG k : " << param.k << "\n";
  param.brk = param_buffer.pop_double();
  console.XDebug()<< "BTriMeshIG r_break: " << param.brk << "\n";
  string buildtype = param_buffer.pop_string();
  console.XDebug()<< "BTriMeshIG build type: " << buildtype.c_str() << "\n";

  // get pointer to mesh
  TriMesh* tmp=NULL;
  if (m_mesh.find(meshname) != m_mesh.end())
  {
    tmp = m_mesh[meshname];
  }
  if(tmp==NULL){
    throw runtime_error("unknown mesh name in TSubLattice<T>::addTriMeshIG:" + meshname);
  }

  // setup new interaction storage
  TriMesh_PIS_EB<T,BTriMeshInteraction>* new_pis=new TriMesh_PIS_EB<T,BTriMeshInteraction>(tmp,m_ppa,param);
  // switch on buildtype, extract buildparam & build
  if(buildtype=="BuildByTag"){
    int tag=param_buffer.pop_int();
    int mask=param_buffer.pop_int();
    new_pis->buildFromPPATagged(tag,mask);
    m_bpis.insert(make_pair(name,new_pis));
  } else if(buildtype=="BuildByGap"){
    double max_gap=param_buffer.pop_double();
    new_pis->buildFromPPAByGap(max_gap);
    m_bpis.insert(make_pair(name,new_pis));
  } else { // unknown type-> throw
    throw runtime_error("unknown build type in TSubLattice<T>::addBondedTriMeshIG:" + buildtype);
  }

  console.XDebug()  << "end TSubLattice<T>::addBondedTriMeshIG()\n";
}

/*!
  Add a 2D mesh. Receive all data from master
*/
template <class T>
void TSubLattice<T>::addMesh2D()
{
  console.XDebug()  << "TSubLattice<T>::addMesh2D()\n";

  // MPI buffer
  CVarMPIBuffer param_buffer(m_comm);
  // data buffers
  vector<MeshNodeData2D> node_recv_buffer;
  vector<MeshEdgeData2D> edge_recv_buffer;

  // receive params
  param_buffer.receiveBroadcast(0);

  // extract name
  string name = param_buffer.pop_string();
  console.XDebug()<< "Mesh2D name: " << name.c_str() << "\n";

  // receive nodes
  m_tml_comm.recv_broadcast_cont_packed(node_recv_buffer,0);
  console.XDebug() << "recvd " << node_recv_buffer.size() << " nodes \n";

  // receive edges
  m_tml_comm.recv_broadcast_cont_packed(edge_recv_buffer,0);
  console.XDebug() << "recvd " << edge_recv_buffer.size() << " edges \n";

  // load mesh into new 2D Mesh
  Mesh2D* new_tm=new Mesh2D();
  new_tm->LoadMesh(node_recv_buffer,edge_recv_buffer);

  m_mesh2d.insert(make_pair(name,new_tm));

  console.XDebug()  << "end TSubLattice<T>::addMesh2D()\n";
}

/*!
  Add a (non-bonded) LineMesh (2D) interaction group
*/

template <class T>
void TSubLattice<T>::addMesh2DIG()
{
console.XDebug()  << "TSubLattice<T>::addMesh2DIG()\n";

  // MPI buffer
  CVarMPIBuffer param_buffer(m_comm);

  // receive params
  param_buffer.receiveBroadcast(0);

  // extract name & type
  string type = param_buffer.pop_string();
  console.XDebug()<< "Mesh2DIG type: " << type.c_str() << "\n";
  string name = param_buffer.pop_string();
  console.XDebug()<< "Mesh2DIG name: " << name.c_str() << "\n";
  string meshname = param_buffer.pop_string();
  console.XDebug()<< "Mesh2DIG mesh name: " << meshname.c_str() << "\n";

  // get pointer to mesh
  Mesh2D* tmp=NULL;
  if (m_mesh2d.find(meshname) != m_mesh2d.end())
  {
    tmp = m_mesh2d[meshname];
  }
  if(tmp==NULL){
    throw runtime_error("unknown mesh name in TSubLattice<T>::addMesh2DIG:" + meshname);
  }
  // switch on type,extract params & construc new Mesh2DIG
  if(type=="Elastic")
  {
    AParallelInteractionStorage* new_pis;
    ETriMeshIP tmi;
    tmi.k=param_buffer.pop_double();
    new_pis = new Mesh2D_PIS_NE<T,EMesh2DInteraction>(tmp,m_ppa,tmi);
    m_dpis.insert(make_pair(name,new_pis));
  } else { // unknown type-> throw
    throw runtime_error("unknown type in TSubLattice<T>::addMesh2DIG:" + type);
  }


  console.XDebug()  << "end TSubLattice<T>::addTriMeshIG()\n";
}

/*!
  add bonded interactions with 2d mesh
  get parameters from master
*/
template <class T>
void TSubLattice<T>::addBondedMesh2DIG()
{
  console.XDebug()  << "TSubLattice<T>::addBondedMesh2DIG()\n";

  // MPI buffer
  CVarMPIBuffer param_buffer(m_comm);

  // receive params
  BMesh2DIP param;
  param_buffer.receiveBroadcast(0);

  // extract name.meshname & params
  string name = param_buffer.pop_string();
  console.XDebug() << "BMesh2DIG name: " << name.c_str() << "\n";
  string meshname = param_buffer.pop_string();
  console.XDebug() << "BMesh2DIG mesh name: " << meshname.c_str() << "\n";
  param.k = param_buffer.pop_double();
  console.XDebug() << "BMesh2DIG k : " << param.k << "\n";
  param.brk = param_buffer.pop_double();
  console.XDebug() << "BMesh2DIG r_break: " << param.brk << "\n";
  string buildtype = param_buffer.pop_string();
  console.XDebug() << "BMesh2DIG build type: " << buildtype.c_str() << "\n";

  // get pointer to mesh
  Mesh2D* tmp=NULL;
  if (m_mesh2d.find(meshname) != m_mesh2d.end())
  {
    tmp = m_mesh2d[meshname];
  }
  if(tmp==NULL){
    throw runtime_error("unknown mesh name in TSubLattice<T>::addBondedMesh2DIG:" + meshname);
  }

  // setup new interaction storage
  Mesh2D_PIS_EB<T,BMesh2DInteraction>* new_pis=new Mesh2D_PIS_EB<T,BMesh2DInteraction>(tmp,m_ppa,param);
  // switch on buildtype, extract buildparam & build
  if(buildtype=="BuildByTag"){
    int tag=param_buffer.pop_int();
    int mask=param_buffer.pop_int();
    new_pis->buildFromPPATagged(tag,mask);
    m_bpis.insert(make_pair(name,new_pis));
  } else if(buildtype=="BuildByGap"){
    double max_gap=param_buffer.pop_double();
    new_pis->buildFromPPAByGap(max_gap);
    m_bpis.insert(make_pair(name,new_pis));
  } else { // unknown type-> throw
    throw runtime_error("unknown build type in TSubLattice<T>::addBondedMesh2DIG:" + buildtype);
  }

  console.XDebug() << "end TSubLattice<T>::addBonded2DMeshIG()\n";
}


/*!
  Add a SingleInteractionGroup to the lattice

  \todo currently gravity is hardwired, change this
*/
template <class T>
void TSubLattice<T>::addSingleIG()
{
  console.XDebug()  << "TSubLattice<T>::addSingleIG()\n";
  CVarMPIBuffer param_buffer(m_comm);

  // get params
  param_buffer.receiveBroadcast(0);

  string type = param_buffer.pop_string();
  console.XDebug()<< "SIG type: " << type.c_str() << "\n";

  // setup InteractionGroup
  if(type=="Gravity"){
    esys::lsm::BodyForceIGP prms = esys::lsm::BodyForceIGP::extract(&param_buffer);

    // add to map
    m_singleParticleInteractions.insert(
      std::pair<string,AInteractionGroup<T>*>(
        prms.Name(),
        new esys::lsm::BodyForceGroup<T>(prms, *m_ppa)
      )
    );
  }
  else if (type=="Buoyancy"){
    esys::lsm::BuoyancyIGP prms = esys::lsm::BuoyancyIGP::extract(&param_buffer);

    // add to map
    m_singleParticleInteractions.insert(
      std::pair<string,AInteractionGroup<T>*>(
        prms.Name(),
        new esys::lsm::BuoyancyForceGroup<T>(prms, *m_ppa)
      )
    );
  }
  else {
    throw std::runtime_error(
      std::string("Trying to setup SIG of unknown type: ")
      +
      type
    );
  }

  console.XDebug()  << "end TSubLattice<T>::addSingleIG()\n";
}


/*!
  Add a DampingGroup to the lattice. Receive the parameters from master.
*/
template <class T>
void TSubLattice<T>::addDamping()
{
  console.XDebug()  << "TSubLattice<T>::addDamping()\n";
  CVarMPIBuffer param_buffer(m_comm);
  // get params
  param_buffer.receiveBroadcast(0);

  string type = param_buffer.pop_string();
  console.XDebug()<< "Damping type: " << type.c_str() << "\n";

  // setup InteractionGroup
  doAddDamping(type,param_buffer);

  console.XDebug()  << "end TSubLattice<T>::addDamping()\n";
}

/*!
  Do the work for adding the damping

  \param type the type of damping
  \param param_buffer the buffer containing the parameters
*/
template <class T>
bool TSubLattice<T>::doAddDamping(const string& type,CVarMPIBuffer& param_buffer)
{
  AParallelInteractionStorage* DG;
  string damping_name;
  bool res;

  if(type=="Damping")
  {
    CDampingIGP *params=extractDampingIGP(&param_buffer);
    DG=new ParallelInteractionStorage_Single<T,CDamping<T> >(m_ppa,*params);
    damping_name="Damping";
    res=true;
  } else if (type=="ABCDamping"){
    ABCDampingIGP *params=extractABCDampingIGP(&param_buffer);
    DG=new ParallelInteractionStorage_Single<T,ABCDamping<T> >(m_ppa,*params);
    damping_name=params->getName();
    res=true;
  } else if (type=="LocalDamping"){
    CLocalDampingIGP *params=extractLocalDampingIGP(&param_buffer);
    DG=new ParallelInteractionStorage_Single<T,CLocalDamping<T> >(m_ppa,*params);
    damping_name=params->getName();
    res=true;
  }else {
    std::stringstream msg;
    msg << "Trying to setup Damping of unknown type: " << type;
    console.Error() << msg.str() << "\n";
    throw std::runtime_error(msg.str());
    res=false;
  }

  // add to map
  if(res) {
    m_damping.insert(make_pair(damping_name,DG));
    m_damping[damping_name]->update();
  }
  return res;
}

/*!
  Add bonded interaction group to the lattice. Receive the parameters from master.
  The bonds are created from the neighbor table.
*/
template <class T>
void TSubLattice<T>::addBondedIG()
{
  console.XDebug()  << "TSubLattice<T>::addBondedIG()\n";
  CVarMPIBuffer param_buffer(m_comm);

  // get params
  CBondedIGP param;
  param_buffer.receiveBroadcast(0);
  param.tag=param_buffer.pop_int();
  string name = param_buffer.pop_string();
  param.k=param_buffer.pop_double();
  param.rbreak=param_buffer.pop_double();
  param.m_scaling=static_cast<bool>(param_buffer.pop_int());

  console.XDebug()
    << "Got BondedIG parameters: " << param.tag
    << " " << name.c_str() << " "
    << param.k << " " << param.rbreak << "\n";
  // setup InteractionGroup
  ParallelInteractionStorage_EB<T,CBondedInteraction> *B_PIS =
    new ParallelInteractionStorage_EB<T,CBondedInteraction>(m_ppa,param);

  // set unbreakable if rbeak<0
  if(param.rbreak<0){
    B_PIS->setUnbreakable(true);
    console.XDebug() << "set bpis unbreakable\"n";
  }

  vector<int> vi(2,-1);
  for(size_t i=0;i<m_temp_conn[param.tag].size();i+=2)
  {
    vi[0] = (m_temp_conn[param.tag][i]);
    vi[1] = (m_temp_conn[param.tag][i+1]);
    B_PIS->tryInsert(vi);
  }

  // add InteractionGroup to map
  m_bpis.insert(make_pair(name,B_PIS));

  console.XDebug()  << "end TSubLattice<T>::addBondedIG()\n";
}

/*!
  Add bonded interaction group to the lattice. Receive the parameters from master.
  The bonds are created from the neighbor table.
*/
template <class T>
void TSubLattice<T>::addCappedBondedIG()
{
  console.XDebug()  << "TSubLattice<T>::addCappedBondedIG()\n";
  CVarMPIBuffer param_buffer(m_comm);

  // get params
  param_buffer.receiveBroadcast(0);
  int tag=param_buffer.pop_int();
  string name = param_buffer.pop_string();
  double k=param_buffer.pop_double();
  double rbreak=param_buffer.pop_double();
  double maxforce=param_buffer.pop_double();

  console.XDebug()
    << "Got CappedBondedIG parameters: " << tag
    << " " << name.c_str() << " "
    << k << " " << rbreak << " " << maxforce << "\n";
  // setup InteractionGroup
  CCappedBondedIGP param;
  param.k=k;
  param.rbreak=rbreak;
  param.tag = tag;
  param.m_force_limit=maxforce;
  ParallelInteractionStorage_EB<T,CCappedBondedInteraction> *B_PIS =
    new ParallelInteractionStorage_EB<T,CCappedBondedInteraction>(m_ppa,param);

  // set unbreakable if rbeak<0
  if(rbreak<0){
    B_PIS->setUnbreakable(true);
    console.XDebug() << "set bpis unbreakable\"n";
  }
  // recv broadcast connection data
  /*console.XDebug()
    << "rank=" << m_tml_comm.rank()
    << "TSubLattice<T>::addCappedBondedIG(): receiving conn_data.\n";

  vector<int> conn_data;
  m_tml_comm.recv_broadcast_cont(conn_data,0);
  console.XDebug()
    << "rank=" << m_tml_comm.rank()
    << "TSubLattice<T>::addBondedIG(): conn_data.size()="
    << conn_data.size() << "\n";
  */
  vector<int> vi(2,-1);
  for(size_t i=0;i<m_temp_conn[tag].size();i+=2)
  {
    vi[0] = (m_temp_conn[tag][i]);
    vi[1] = (m_temp_conn[tag][i+1]);
    B_PIS->tryInsert(vi);
  }

  // add InteractionGroup to map
  m_bpis.insert(make_pair(name,B_PIS));

  console.XDebug()  << "end TSubLattice<T>::addCappedBondedIG()\n";
}

template <class T>
void TSubLattice<T>::addRotBondedIG()
{
  console.Error()  << "TSubLattice<T>::addRotBondedIG() => trying to add rotational bonded IG to nonrotational model\n";
}

template <class T>
void TSubLattice<T>::addRotThermBondedIG()
{
  console.Error()  << "TSubLattice<T>::addRotThermBondedIG() => trying to add rotational thermal bonded IG to nonrotational model\n";
}

template <class T>
void TSubLattice<T>::addBrittleBeamSCIG()
{
    console.Error()  << "TSubLattice<T>::addBrittleBeamSCIG() => trying to add rotational bonded IG to nonrotational model\n";
}

template <class T>
void TSubLattice<T>::addBrittleBeamDZCIG()
{
    console.Error()  << "TSubLattice<T>::addBrittleBeamDZCIG() => trying to add rotational bonded IG to nonrotational model\n";
}

/*!
  Add short bonded interaction group to the lattice. Receive the parameters from master.
  The bonds are created from the neighbor table.
*/
template <class T>
void TSubLattice<T>::addShortBondedIG()
{
  console.XDebug()  << "TSubLattice<T>::addShortBondedIG()\n";
  CVarMPIBuffer param_buffer(m_comm);

  // get params
  param_buffer.receiveBroadcast(0);
  int tag=param_buffer.pop_int();
  string name = param_buffer.pop_string();
  double k=param_buffer.pop_double();
  double rbreak=param_buffer.pop_double();

  console.XDebug()
    << "Got ShortBondedIG parameters: " << tag
    << " " << name.c_str() << " "
    << k << " " << rbreak << "\n";
  // setup InteractionGroup
  CBondedIGP param;
  param.k=k;
  param.rbreak=rbreak;
  param.tag = tag;
  ParallelInteractionStorage_EB<T,CShortBondedInteraction> *B_PIS =
    new ParallelInteractionStorage_EB<T,CShortBondedInteraction>(m_ppa,param);

  // recv broadcast connection data
  /*console.XDebug()
    << "rank=" << m_tml_comm.rank()
    << "TSubLattice<T>::addShortBondedIG(): receiving conn_data.\n";

  vector<int> conn_data;
  m_tml_comm.recv_broadcast_cont(conn_data,0);
  console.XDebug()
    << "rank=" << m_tml_comm.rank()
    << "TSubLattice<T>::addShortBondedIG(): conn_data.size()="
    << conn_data.size() << "\n";*/

  vector<int> vi(2,-1);
  for(size_t i=0;i<m_temp_conn[param.tag].size();i+=2)
  {
    vi[0] = (m_temp_conn[param.tag][i]);
    vi[1] = (m_temp_conn[param.tag][i+1]);
    B_PIS->tryInsert(vi);
  }

  // add InteractionGroup to map
  m_bpis.insert(make_pair(name,B_PIS));

  console.XDebug()  << "end TSubLattice<T>::addShortBondedIG()\n";
}

/*!
  Set excluding IG. The names of 2 interaction groups are received from the master.
  The 1st is set to be excluding the 2nd. If one of the named IGs does not exist,
  the operation is ignored and a warning emitted.
*/
template <class T>
void TSubLattice<T>::setExIG()
{
  console.XDebug()  << "TSubLattice<T>::addExIG()\n";
  CVarMPIBuffer pbuffer(m_comm);

  // get params
  pbuffer.receiveBroadcast(0);
  string s1 = pbuffer.pop_string();
  string s2 = pbuffer.pop_string();

  //console.XDebug()<< s1.c_str()  << "  " <<  s2.c_str() << "\n";
  // for first (excluding) IG, look in bonded _and_ in dynamic
  map<string,AParallelInteractionStorage*>::iterator excluding_ig=m_bpis.find(s1);
  if(excluding_ig==m_bpis.end()){
		excluding_ig=m_dpis.find(s2);
  }
  map<string,AParallelInteractionStorage*>::iterator dynamic_ig=m_dpis.find(s2);
  if((excluding_ig!=m_dpis.end())&&(dynamic_ig!=m_dpis.end()))
  {
    // map iterators point to [key,value] pairs, therefore it->second
    // is the pointer to the PIS here
    dynamic_ig->second->addExIG(excluding_ig->second);
  }
  else
  {
    console.Error() << "TSubLattice<T>::setExIG() - nonexisting interaction group \n";
  }

  console.XDebug()  << "end TSubLattice<T>::addExIG()\n";
}

/*!
  Remove interaction group. The name of the group is received from the master.

  \warning Doesn't deal yet with removal of dependent items, i.e. savers
*/
template <class T>
void TSubLattice<T>::removeIG()
{
  console.XDebug()  << "TSubLattice<T>::removeIG()\n";
  CVarMPIBuffer pbuffer(m_comm);
  bool found=false;

  // get params
  pbuffer.receiveBroadcast(0);
  string igname = pbuffer.pop_string();
  // look for name in map of non-bondend particle pair interactions
  map<string,AParallelInteractionStorage*>::iterator iter=m_dpis.find(igname);
  // if found, remove interaction
  if(iter!=m_dpis.end()){
    found=true;
    delete m_dpis[igname];
    m_dpis.erase(iter);
  } else {
    // if not found, look in unbonded wall interactions
    typename map<string,AWallInteractionGroup<T>*>::iterator it2=m_WIG.find(igname);
    if(it2!=m_WIG.end()){
      found=true;
      delete m_WIG[igname];
      m_WIG.erase(it2);
    }
  }

  if(!found) {
    console.Error() << "TSubLattice<T>::removeIG() - nonexisting interaction group - ignore removal\n";
  }
}


/*!
  Exchange positions of shared particles with the other sublattices.
*/
template <class T>
void TSubLattice<T>::exchangePos()
{
  console.XDebug() << "TSubLattice<T>::exchangePos() \n" ;

  m_ppa->exchange(&T::getExchangeValues,&T::setExchangeValues);

  console.XDebug() << "end TSubLattice<T>::exchangePos()\n" ;
}

/*!
  Reset the forces on all particles and walls to 0
*/
template <class T>
void TSubLattice<T>::zeroForces()
{
  console.XDebug() << "TSubLattice<T>::zeroForces()\n";

  // particles
  m_ppa->forAllParticles(&T::zeroForce);
  // trimeshes
  for(map<string,TriMesh*>::iterator iter=m_mesh.begin();
      iter!=m_mesh.end();
      iter++){
    (iter->second)->zeroForces();
  }

  // 2d meshes
  for(map<string,Mesh2D*>::iterator iter=m_mesh2d.begin();
      iter!=m_mesh2d.end();
      iter++){
    (iter->second)->zeroForces();
  }
  // walls
  for(typename map<string,CWall*>::iterator iter=m_walls.begin();
      iter!=m_walls.end();
      iter++)
  {
    (iter->second)->zeroForce();
  }
  // sphere bodies
  for(typename map<string,CSphereBody*>::iterator iter=m_spheres.begin();
      iter!=m_spheres.end();
      iter++)
  {
    (iter->second)->zeroForce();
  }
  console.XDebug() << "end TSubLattice<T>::zeroForces() \n";
}

/*!
  Calculate the forces for all interactions contained in the sublattice.
  Interactions contained in more than one sublattice are calculated in
  each of them. Slightly more computation, but less communication, i.e.
  saves one syncronisation point (exchange of forces).
*/
template <class T>
void TSubLattice<T>::calcForces()
{
  console.XDebug() << "TSubLattice<T>::calcForces() \n";


  // the walls
  for(typename map<string,AWallInteractionGroup<T>*>::iterator it=m_WIG.begin();it!=m_WIG.end();it++)
  {
    (it->second)->calcForces();
  }
  // the sphere bodies
  for(typename map<string,ASphereBodyInteractionGroup<T>*>::iterator it=m_SIG.begin();it!=m_SIG.end();it++)
  {
    (it->second)->calcForces();
  }
  // single particle IGs
  for(
    typename NameIGroupMap::iterator siter=this->m_singleParticleInteractions.begin();
    siter != m_singleParticleInteractions.end();
    siter++
  )
  {
    (siter->second)->calcForces();
  }
  // dynamically created IGs
  for(typename map<string,AParallelInteractionStorage*>::iterator iter=m_dpis.begin();iter!=m_dpis.end();iter++)
  {
    (iter->second)->calcForces();
  }
  // bonded IGs
  for(typename map<string,AParallelInteractionStorage*>::iterator iter=m_bpis.begin();iter!=m_bpis.end();iter++)
  {
    (iter->second)->calcForces();
  }
  // last, but not least - damping
  for(typename map<string,AParallelInteractionStorage*>::iterator iter=m_damping.begin();iter!=m_damping.end();iter++)
  {
    (iter->second)->calcForces();
  }


  if(fluidExist()){m_fluidinteraction->calcForces();}; //fluid contents
  

  console.XDebug() << "end TSubLattice<T>::calcForces() \n";
}

/*!
  Calculate the forces for all interactions contained in the sublattice.
  Interactions contained in more than one sublattice are calculated in
  each of them. Slightly more computation, but less communication, i.e.
  saves one syncronisation point (exchange of forces).
*/
template <class T>
void TSubLattice<T>::setTimeStepSize(double dt)
{
  m_dt = dt;
  console.XDebug() << "TSubLattice<T>::setTimeStepSize() \n";

  // the walls
  for(typename map<string,AWallInteractionGroup<T>*>::iterator it=m_WIG.begin();it!=m_WIG.end();it++)
  {
    (it->second)->setTimeStepSize(dt);
  }
  // the sphere bodies
  for(typename map<string,ASphereBodyInteractionGroup<T>*>::iterator it=m_SIG.begin();it!=m_SIG.end();it++)
  {
    (it->second)->setTimeStepSize(dt);
  }
  // single particle IGs
  for(
    typename NameIGroupMap::iterator siter=this->m_singleParticleInteractions.begin();
    siter != m_singleParticleInteractions.end();
    siter++
  )
  {
    (siter->second)->setTimeStepSize(dt);
  }
  // dynamically created IGs
  for(typename map<string,AParallelInteractionStorage*>::iterator iter=m_dpis.begin();iter!=m_dpis.end();iter++)
  {
    (iter->second)->setTimeStepSize(dt);
  }
  // bonded IGs
  for(typename map<string,AParallelInteractionStorage*>::iterator iter=m_bpis.begin();iter!=m_bpis.end();iter++)
  {
    (iter->second)->setTimeStepSize(dt);
  }
  // last, but not least - damping
  for(typename map<string,AParallelInteractionStorage*>::iterator iter=m_damping.begin();iter!=m_damping.end();iter++)
  {
    (iter->second)->setTimeStepSize(dt);
  }

  console.XDebug() << "end TSubLattice<T>::setTimeStepSize() \n";
}

/*!
  Do the time integration for the particles owned by the SubLattice by 1st order method

  \param dt the time step
 */
template <class T>
void TSubLattice<T>::integrate(double dt)
{
  console.XDebug() << "TSubLattice<T>::integrate \n";
  m_ppa->forAllParticles(&T::integrate,dt);
  m_ppa->forAllParticles(&T::rescale) ;
  console.XDebug() << "end TSubLattice<T>::integrate \n";
}

/*!
  Do one step, i.e. calculate forces, velocities, positions
*/
template <class T>
void TSubLattice<T>::oneStep()
{
  zeroForces();
  m_pTimers->start("TSubLatticecalcForces");
  calcForces();
  m_pTimers->stop("TSubLatticecalcForces");
  integrate(m_dt);

  if (this->getParticleType() == "RotTherm")
  {
    this->oneStepTherm();
  }
}

/*!
  Do one step, i.e. calculate forces, velocities, positions
*/
template <class T>
void TSubLattice<T>::oneStepTherm()
{
  zeroHeat();            // ???? combine?
  calcHeatFrict();
  calcHeatTrans();
  integrateTherm(m_dt);
  thermExpansion() ;
}

/*!
  Do the time integration for the particles owned by the SubLattice for temprature

  \param dt the time step
 */
template <class T>
void TSubLattice<T>::integrateTherm(double dt)
{
  console.XDebug() << "TSubLattice<T>::integrateTherm \n";
  m_ppa->forAllParticles(&T::integrateTherm,dt);
//  m_ppa->forAllParticles(&T::rescale) ;
  console.XDebug() << "end TSubLattice<T>::integrateTherm \n";
}

template <class T>
void TSubLattice<T>::thermExpansion()
{
  console.XDebug() << "TSubLattice<T>::thermExpansion() \n";
  m_ppa->forAllParticles(&T::thermExpansion);
//  m_ppa->forAllParticles(&T::rescale) ;
  console.XDebug() << "end TSubLattice<T>::thermExpansion() \n";
}

/*!
  Reset the HeatSources on all particles and walls to 0
*/
template <class T>
void TSubLattice<T>::zeroHeat()
{
  console.XDebug() << "TSubLattice<T>::zeroHeat()\n";

  // particles
  m_ppa->forAllParticles(&T::zeroHeat);

/*
// walls
  for(typename vector<AWallInteractionGroup<T>*>::iterator iter=m_WIG.begin();iter!=m_WIG.end();iter++)
  {
    (*iter)->zeroForce();
  }
*/
  console.XDebug() << "end TSubLattice<T>::zeroHeat() \n";
}

/*!
  Calculate the Heat Sources for all interactions contained in the sublattice.
  Interactions contained in more than one sublattice are calculated in
  each of them. Slightly more computation, but less communication, i.e.
  saves one syncronisation point (exchange of forces).
*/
template <class T>
void TSubLattice<T>::calcHeatFrict()
{
  console.XDebug() << "TSubLattice<T>::calcHeatFrict() \n";

  // dynamically created IGs
  for(
    typename map<string,AParallelInteractionStorage*>::iterator iter=m_dpis.begin();
    iter!=m_dpis.end();
    iter++
  )
  {
    (iter->second)->calcHeatFrict();
  }

  console.XDebug() << "end TSubLattice<T>::calcHeatFrict() \n";
}

template <class T>
void TSubLattice<T>::calcHeatTrans()
{
  console.XDebug() << "TSubLattice<T>::calcHeatTrans() \n";


  // dynamically created IGs
  for(
    typename map<string,AParallelInteractionStorage*>::iterator iter=m_dpis.begin();
    iter!=m_dpis.end();
    iter++
  )
  {
    (iter->second)->calcHeatTrans();
  }
  // bonded IGs
  for(
    typename map<string,AParallelInteractionStorage*>::iterator iter=m_bpis.begin();
    iter!=m_bpis.end();
    iter++
  )
  {
    (iter->second)->calcHeatTrans();
  }

  console.XDebug() << "end TSubLattice<T>::calcHeatTrans() \n";
}

/*!
  rebuild particle array
*/
template <class T>
void TSubLattice<T>::rebuildParticleArray()
{
  m_ppa->rebuild();
}

/*!
  rebuild interactions
*/
template <class T>
void TSubLattice<T>::rebuildInteractions()
{
  CMPIBarrier barrier(m_worker_comm);
  m_pTimers->start("RebuildInteractions");
  m_pTimers->resume("NeighbourSearch");
  for(typename map<string,AParallelInteractionStorage*>::iterator iter=m_bpis.begin();
          iter!=m_bpis.end();
          iter++)
  {
    console.Debug() << "exchg & rebuild BPIS " << iter->first << " at node " << m_rank << "\n";
    (iter->second)->exchange();
    (iter->second)->rebuild();
  }

  for(typename map<string,AParallelInteractionStorage*>::iterator iter=m_dpis.begin();
      iter!=m_dpis.end();
      iter++)
  {
    console.Debug() << "exchg & rebuild DPIS " << iter->first << " at node " << m_rank << "\n";
    (iter->second)->exchange();
    m_pTimers->pause("RebuildInteractions");
    m_pTimers->pause("NeighbourSearch");
    barrier.wait("dpis::exchange");
    m_pTimers->resume("RebuildInteractions");
    m_pTimers->resume("NeighbourSearch");
    (iter->second)->rebuild();
  }
  resetDisplacements();
  m_pTimers->stop("RebuildInteractions");
}

/*!
  Search neighbors - rebuild particle array and interactions
*/
template <class T>
void TSubLattice<T>::searchNeighbors()
{
  console.Debug() << "CSubLattice<T>::searchNeighbors()\n";
  CMPIBarrier barrier(getWorkerComm());
  m_pTimers->start("NeighbourSearch");
  m_pTimers->start("RebuildParticleArray");
  rebuildParticleArray();
  m_pTimers->stop("RebuildParticleArray");
  m_pTimers->pause("NeighbourSearch");
  barrier.wait("PPA rebuild");
  rebuildInteractions();
  m_pTimers->stop("NeighbourSearch");
  console.Debug() << "end CSubLattice<T>::searchNeighbors()\n";
}

/*!
  Update the interaction groups from an existing Neighbortable

  \todo check for rebuild Neighbortable
*/
template <class T>
void TSubLattice<T>::updateInteractions()
{
  console.Debug() << "updateInteractions() \n";
  console.Debug() << "m_ppa->getTimeStamp() " << m_ppa->getTimeStamp() << " m_last_ns " << m_last_ns << "\n";
  bool need_update=false;

  m_pTimers->start("UpdateBondedInteractions");
  for(typename map<string,AParallelInteractionStorage*>::iterator iter=m_bpis.begin();
          iter!=m_bpis.end();
          iter++)
  {
    bool n=(iter->second)->update();
    need_update=need_update || n;
  }
  m_pTimers->stop("UpdateBondedInteractions");
  if((m_ppa->getTimeStamp() > m_last_ns) || need_update)
  {
    m_pTimers->start("UpdateDynamicInteractions");
    for(typename map<string,AParallelInteractionStorage*>::iterator iter=m_dpis.begin();
          iter!=m_dpis.end();
          iter++)
      {
	bool n=(iter->second)->update();
	need_update=need_update || n;
      }
    m_pTimers->stop("UpdateDynamicInteractions");
    for(typename map<string,AWallInteractionGroup<T>*>::iterator it=m_WIG.begin();
        it!=m_WIG.end();
        it++)
      {
	(it->second)->Update(m_ppa);
      }
    for(typename map<string,ASphereBodyInteractionGroup<T>*>::iterator it=m_SIG.begin();
        it!=m_SIG.end();
        it++)
      {
	(it->second)->Update(m_ppa);
      }
    for(typename map<string,AParallelInteractionStorage*>::iterator iter=m_damping.begin();
         iter!=m_damping.end();
         iter++){
      (iter->second)->update();
    }
    
    if(fluidExist()){m_fluidinteraction->update();};//fluid contents

    m_last_ns=m_ppa->getTimeStamp();
  }

  console.Debug() << "end TSubLattice<T>::updateInteractions()\n";
}

/*!
  check if any of the owned particles has moved futher than the search range for the neighbor table
  5 Flops(1 dotproduct)/particle looked at
*/
template <class T>
void TSubLattice<T>::checkNeighbors()
{
  console.Debug() << "TSubLattice<T>::checkNeighbors()\n";
  CMPISGBufferLeaf buffer(m_comm,0,8);
  double mdsqr=0; // squared max. displacement
  double alpha=0.5*m_alpha;
  double srsqr=alpha*alpha; // squared search range
  vector<Vec3> displ; // displacements
  int res; // result 0/1

  // --- particles ---
  // get displacement data
  m_ppa->forAllParticlesGet(displ,&T::getDisplacement);

  // find maximum particle displacement
  vector<Vec3>::iterator it=displ.begin();
  while((it!=displ.end())&&(mdsqr<srsqr))
  {
    double sqdisp=(*it)*(*it);
    mdsqr = ((mdsqr < sqdisp) ? sqdisp : mdsqr);
    it++;
    //console.XDebug() << "sq. disp: " << sqdisp << "\n";
  }
  console.XDebug() << "max squared displacement " << mdsqr << " alpha^2 = " << srsqr << "\n";
  if (mdsqr>srsqr){
    res=1;
  } else {
    res=0;
  }

  // --- mesh ---
  // only needed if res==0
  if(res==0){
    for(map<string,TriMesh*>::iterator iter=m_mesh.begin();
	iter!=m_mesh.end();
	iter++){
      console.XDebug() << "checking mesh " << iter->first << "\n";
      if(iter->second->hasMovedBy(alpha)){
		res=1;
	console.XDebug() << "mesh has moved too far\n";
	}
    }
  }
  buffer.append(res);
  buffer.send();
  console.Debug() << "end TSubLattice<T>::checkNeighbors()\n";
}


/*!
  Reset the displacement of all particles & meshes
*/
template <class T>
void TSubLattice<T>::resetDisplacements()
{
  console.Debug() << "slave " << m_rank << " resetDisplacements()\n";
  m_ppa->forAllParticles(&T::resetDisplacement);
  for(map<string,TriMesh*>::iterator iter=m_mesh.begin();
      iter!=m_mesh.end();
      iter++){
    iter->second->resetCurrentDisplacement();
  }
  console.Debug() << "slave " << m_rank << " end resetDisplacements()\n";
}

/*!
  Change the radius of all particles with a given tag by a specified amount.
  Parameters (tag,deltaR) are received from master.
*/
template <class T>
void TSubLattice<T>::changeRadiusBy()
{
  console.Debug() << "TSubLattice<T>::changeRadiusBy()\n";
  CVarMPIBuffer buffer(m_comm);

  buffer.receiveBroadcast(0); // get data from master
  int tag=buffer.pop_int();
  double deltaR=buffer.pop_double();
  m_ppa->forParticleTag(tag,(void (T::*)(double))(&T::changeRadiusBy),deltaR);
  console.Debug() << "end TSubLattice<T>::changeRadiusBy()\n";
}

/*!
  Move all particles with a given tag to a given position.
  Parameters (tag,posn) are received from master.
*/
template <class T>
void TSubLattice<T>::moveParticleTo()
{
  console.Debug() << "TSubLattice<T>::moveParticleTo()\n";
  CVarMPIBuffer buffer(m_comm);

  buffer.receiveBroadcast(0); // get data from master
  int tag=buffer.pop_int();
  Vec3 mv=buffer.pop_vector();
  m_ppa->forParticleTag(tag,(void (T::*)(Vec3))(&T::moveToRel),mv);
  console.Debug() << "end TSubLattice<T>::moveParticleTo()\n";
}

/*!
  Move all particles with a given tag by a specified displacement.
  Parameters (tag,displacement) are received from master.
*/
template <class T>
void TSubLattice<T>::moveTaggedParticlesBy()
{
  console.Debug() << "TSubLattice<T>::moveTaggedParticlesBy()\n";
  CVarMPIBuffer buffer(m_comm);

  buffer.receiveBroadcast(0); // get data from master
  const int tag = buffer.pop_int();
  const Vec3 dx = buffer.pop_vector();
  m_ppa->forParticleTag(tag, (void (T::*)(Vec3))(&T::moveBy),dx);
  console.Debug() << "end TSubLattice<T>::moveTaggedParticlesBy()\n";
}


template <class T>
void TSubLattice<T>::moveSingleParticleTo(int particleId, const Vec3 &posn)
{
  m_ppa->forParticle(particleId, (void (T::*)(Vec3))(&T::moveTo), posn);
}

/*!
  Move mesh (tri or 2d) node by a given amount.
  Parameters (name,id,displacement) are received from master.
*/
template <class T>
void TSubLattice<T>::moveSingleNode()
{
  console.Debug() << "TSubLattice<T>::moveSingleNode()\n";
  CVarMPIBuffer buffer(m_comm);

  buffer.receiveBroadcast(0); // get data from master
  string name=string(buffer.pop_string());
  int id=buffer.pop_int();
  Vec3 disp=buffer.pop_vector();

  console.XDebug() << "name :" << name << " id : " << id << " disp " << disp << "\n";

  map<string,TriMesh*>::iterator tm=m_mesh.find(name);
  if (tm!=m_mesh.end()){
    (tm->second)->moveNode(id,disp);
  } else {
    map<string,Mesh2D*>::iterator m2d=m_mesh2d.find(name);
    if(m2d!=m_mesh2d.end()){
      (m2d->second)->moveNode(id,disp);
    }
  }
  console.Debug() << "end TSubLattice<T>::moveSingleNode()\n";
}

/*!
  Move tagged trimesh nodes by a given amount.
  Parameters (name,tag,displacement) are received from master.
*/
template <class T>
void TSubLattice<T>::moveTaggedNodes()
{
  console.Error() << "TSubLattice<T>::moveTaggedNodes() NOT IMPLEMENTED\n";
  throw
    std::runtime_error(
      "TSubLattice<T>::moveTaggedNodes() NOT IMPLEMENTED\n"
    );
#if 0
  CVarMPIBuffer buffer(m_comm);

  buffer.receiveBroadcast(0); // get data from master
  string name=string(buffer.pop_string());
  int tag=buffer.pop_int();
  Vec3 disp=buffer.pop_vector();
#endif
}

/*!
  translate mesh by given amount

  \param meshName the name of the mesh to be moved
  \param translation the translation vector
*/
template <class T>
void TSubLattice<T>::translateMeshBy(const std::string &meshName, const Vec3 &translation)
{
  map<string,TriMesh*>::iterator tm=m_mesh.find(meshName);
  if (tm != m_mesh.end()){
    (tm->second)->translateBy(translation);
  } 
  map<string,Mesh2D*>::iterator m2d=m_mesh2d.find(meshName);
  if(m2d!=m_mesh2d.end()){
    (m2d->second)->translateBy(translation);
  }
}

/*!
    rotate mesh around axis by given amount

    \param meshName the name of the mesh to be moved
    \param origin a point on the axis
    \param axis the orientation of the rotation axis
    \param angle the rotation angle in degrees
*/
template <class T>
void TSubLattice<T>::rotateMeshBy(const std::string &meshName, const Vec3& origin, const Vec3 &axis, double angle)
{
    console.Debug() << "rotating mesh " << meshName << "by (" << axis << " ), " << angle << "\n";  
    map<string,TriMesh*>::iterator tm=m_mesh.find(meshName); // only for 3D meshes!
    console.XDebug() << "past mesh.find\n";
    if (tm != m_mesh.end()){
        console.XDebug() << "found mesh, calling rotateBy\n";
        (tm->second)->rotateBy(origin, axis,angle);
    }  else {
        console.Error() << "trying to rotate non-existing mesh " << meshName << "\n";
    }
}

  

template <class T>
std::pair<double, int> TSubLattice<T>::findParticleNearestTo(const Vec3 &pt)
{
  console.Debug() << "TSubLattice<T>::findParticleNearestTo: enter\n";
  const T *pClosest = NULL;
  double minDistSqrd = std::numeric_limits<double>::max();

  typename ParticleArray::ParticleIterator it =
      m_ppa->getInnerParticleIterator();
  while (it.hasNext())
  {
    const T &p = it.next();
    const double distSqrd = (pt - p.getPos()).norm2();
    if (distSqrd < minDistSqrd)
    {
      minDistSqrd = distSqrd;
      pClosest = &p;
    }
  }
  console.Debug() << "TSubLattice<T>::findParticleNearestTo: exit\n";
  return
    (
      (pClosest != NULL)
      ?
      std::make_pair(sqrt(minDistSqrd), pClosest->getID())
      :
      std::make_pair(std::numeric_limits<double>::max(), -1)
    );
}

/*!
  \todo DO DOCUMENTATION
*/
template <class T>
std::pair<int, Vec3> TSubLattice<T>::getParticlePosn(int particleId)
{
  const T *particle = NULL;
  typename ParticleArray::ParticleIterator it =
      m_ppa->getInnerParticleIterator();
  while (it.hasNext())
  {
    const T &p = it.next();
    if (p.getID() == particleId)
    {
      particle = &p;
    }
  }
  if (particle != NULL)
  {
    return std::make_pair(particleId, particle->getPos());
  }
  return std::make_pair(-1,Vec3::ZERO);
}

/*!
  \todo DO DOCUMENTATION
*/
template <class T>
void TSubLattice<T>::getParticleData(const IdVector &particleIdVector)
{
  console.Debug()
    << "TSubLattice<T>::getParticleData: enter\n";
  typedef std::set<int> IdSet;
  typedef std::vector<T> ParticleVector;

  ParticleVector particleVector;
  typename ParticleArray::ParticleIterator it =
      m_ppa->getInnerParticleIterator();
  if (particleIdVector.size() > 0)
  {
    IdSet idSet(particleIdVector.begin(), particleIdVector.end());
    console.Debug()
      << "TSubLattice<T>::getParticleData: iterating over particles\n";
    while (it.hasNext())
    {
      const T &p = it.next();
      if (idSet.find(p.getID()) != idSet.end())
      {
        particleVector.push_back(p);
      }
    }
  }
  else
  {
    m_ppa->getAllInnerParticles(particleVector);
  }
  console.Debug()
    << "TSubLattice<T>::getParticleData:"
    << " sending particle data of size " << particleVector.size() << "\n";
  m_tml_comm.send_gather_packed(particleVector, 0);
  console.Debug()
    << "TSubLattice<T>::getParticleData: exit\n";
}

/*!
  Tag particle closest to given position. Params received from master.
*/
template <class T>
void TSubLattice<T>::tagParticleNearestTo()
{
  console.Debug() << "TSubLattice<T>::tagParticleNearestTo()\n";
  CVarMPIBuffer buffer(m_comm);

  buffer.receiveBroadcast(0); // get data from master
  int tag=buffer.pop_int();
  int mask=buffer.pop_int();
  Vec3 pos=buffer.pop_vector();

  // warning - this is ugly
  T* part_ptr=m_ppa->getParticlePtrByPosition(pos);
  if(part_ptr!=NULL){
    int old_tag=part_ptr->getTag();
    int new_tag=(old_tag & (~mask)) | (tag & mask);
    part_ptr->setTag(new_tag);

    cout << "pos, realpos: " << pos << " " << part_ptr->getPos() << " old tag, new tag " << old_tag << " " << part_ptr->getTag() << endl;
  }
  console.Debug() << "end TSubLattice<T>::tagParticleNearestTo()\n";
}

/*!
  Make tagged particles non-dynamic  i.e. don't update velocity (rot+lin).
  Parameters (tag) are received from master.
*/
template <class T>
void TSubLattice<T>::setParticleNonDynamic()
{
  console.Debug() << "TSubLattice<T>::setParticleNonDynamic()\n";
  CVarMPIBuffer buffer(m_comm);

  buffer.receiveBroadcast(0); // get data from master
  int tag=buffer.pop_int();
  m_ppa->forParticleTag(tag,(void (T::*)())(&T::setNonDynamic));
  console.Debug() << "end TSubLattice<T>::setParticleNonDynamic()\n";
}

/*!
  Make tagged particles non-rotational, i.e. don't update rotational velocity & position
  Parameters (tag) are received from master.
*/
template <class T>
void TSubLattice<T>::setParticleNonRot()
{
  console.Debug() << "TSubLattice<T>::setParticleNonRot()\n";
  CVarMPIBuffer buffer(m_comm);

  buffer.receiveBroadcast(0); // get data from master
  int tag=buffer.pop_int();
  m_ppa->forParticleTag(tag,(void (T::*)())(&T::setNonDynamicRot));
  console.Debug() << "end TSubLattice<T>::setParticleNonRot()\n";
}

/*!
  Make tagged particles lieanr non-dynamic, i.e. don't update linear velocity
  Parameters (tag) are received from master.
*/
template <class T>
void TSubLattice<T>::setParticleNonTrans()
{
  console.Debug() << "TSubLattice<T>::setParticleNonTrans()\n";
  CVarMPIBuffer buffer(m_comm);

  buffer.receiveBroadcast(0); // get data from master
  int tag=buffer.pop_int();
  m_ppa->forParticleTag(tag,(void (T::*)())(&T::setNonDynamicLinear));
  console.Debug() << "end TSubLattice<T>::setParticleNonTrans()\n";
}

/*!
   Set the velocity of a tagged group of particles. Parameters are received from master.
*/
template <class T>
void TSubLattice<T>::setTaggedParticleVel()
{
  console.Debug() << "TSubLattice<T>::setTaggedParticleVel()\n";
  CVarMPIBuffer buffer(m_comm);

  buffer.receiveBroadcast(0); // get data from master
  int tag=buffer.pop_int();
  Vec3 v=buffer.pop_vector();
  m_ppa->forParticleTag(tag,(void (T::*)(Vec3))(&T::setVel),v);
  console.XDebug() << "end TSubLattice<T>::setTaggedParticleVel()\n";
}

/*!
  Move wall by a given vector. Parameters are received from master.
*/
template <class T>
void TSubLattice<T>::moveWallBy()
{
  console.XDebug() << "TSubLattice<T>::moveWallBy()\n";
  CVarMPIBuffer buffer(m_comm);

  buffer.receiveBroadcast(0); // get data from master
  string wname=buffer.pop_string();
  Vec3 mv=buffer.pop_vector();
  typename map<string,CWall*>::iterator iter=m_walls.find(wname);
  if(iter!=m_walls.end())
  {
    (iter->second)->moveBy(mv);
  }
}

/*!
  Move sphere body by a given vector. Parameters are received from master.
*/
template <class T>
void TSubLattice<T>::moveSphereBodyBy()
{
  console.XDebug() << "TSubLattice<T>::moveSphereBodyBy()\n";
  CVarMPIBuffer buffer(m_comm);

  buffer.receiveBroadcast(0); // get data from master
  string wname=buffer.pop_string();
  Vec3 mv=buffer.pop_vector();
  typename map<string,CSphereBody*>::iterator iter=m_spheres.find(wname);
  if(iter!=m_spheres.end())
  {
    (iter->second)->moveBy(mv);
  }
}

/*!
  Change wall orientation by a given vector. Parameters are received from master.
*/
template <class T>
void TSubLattice<T>::setWallNormal()
{
  console.XDebug() << "TSubLattice<T>::setWallNormal()\n";
  CVarMPIBuffer buffer(m_comm);

  buffer.receiveBroadcast(0); // get data from master
  string wname=buffer.pop_string();
  Vec3 wn=buffer.pop_vector();
  typename map<string,CWall*>::iterator iter=m_walls.find(wname);
  if(iter!=m_walls.end())
  {
    (iter->second)->setNormal(wn);
  }
}

/*!
  Apply given force to a wall. Parameters are received from master.
*/
template <class T>
void TSubLattice<T>::applyForceToWall()
{
  console.XDebug() << "TSubLattice<T>::applyForceToWall()\n";
  CVarMPIBuffer buffer(m_comm);

  buffer.receiveBroadcast(0); // get data from master
  string wname=buffer.pop_string();
  Vec3 f=buffer.pop_vector();
  typename map<string,AWallInteractionGroup<T>*>::iterator iter=m_WIG.find(wname);
  if(iter!=m_WIG.end())
  {
    (iter->second)->applyForce(f);
  }
}

/*!
  Set the velocity of a wall. Only for visc. walls, i.e. m_vel is set,
  but position is not affected. Parameters are received from master.
*/
template <class T>
void TSubLattice<T>::setVelocityOfWall()
{
  console.XDebug() << "TSubLattice<T>::setVelocityOfWall()\n";
  CVarMPIBuffer buffer(m_comm);

  buffer.receiveBroadcast(0); // get data from master
  string wname=buffer.pop_string();
  Vec3 v=buffer.pop_vector();
  typename map<string,AWallInteractionGroup<T>*>::iterator iter=m_WIG.find(wname);
  if(iter!=m_WIG.end())
  {
    (iter->second)->setVelocity(v);
  }
}

/*!
  Set the velocity of a particle. Parameters are received from master.
*/
template <class T>
void TSubLattice<T>::setParticleVelocity()
{
  console.Debug() << "TSubLattice<T>::setParticleVelocity()\n";
  CVarMPIBuffer buffer(m_comm);

  buffer.receiveBroadcast(0); // get data from master
  int id=buffer.pop_int();
  Vec3 mv=buffer.pop_vector();
  m_ppa->forParticle(id,(void (T::*)(Vec3))(&T::setVel),mv);
  console.XDebug() << "end TSubLattice<T>::setParticleVelocity()\n";
}

/*!
  Set the density of a tagged group of particles. Parameters are received from master.
*/
template <class T>
void TSubLattice<T>::setParticleDensity()
{
  console.Debug() << "TSubLattice<T>::setParticleDensity()\n";
  CVarMPIBuffer buffer(m_comm);

  buffer.receiveBroadcast(0); // get data from master
  int tag=buffer.pop_int();
  int mask=buffer.pop_int();
  double rho=buffer.pop_double();
  m_ppa->forParticleTagMask(tag,mask,(void (T::*)(double))(&T::setDensity),rho);
  console.XDebug() << "end TSubLattice<T>::setParticleVelocity()\n";
}

/*!
  reset the orientation of a tagged group of particles. Parameters are received from master.
*/
template <class T>
void TSubLattice<T>::resetParticleRotation()
{
	console.Debug() << "TSubLattice<T>::resetParticleRotation()\n";
    CVarMPIBuffer buffer(m_comm);

    buffer.receiveBroadcast(0); // get data from master
    int tag=buffer.pop_int();
    int mask=buffer.pop_int();
    double rho=buffer.pop_double();
    m_ppa->forParticleTagMask(tag,mask,(void (T::*)())(&T::resetRotation));
    console.XDebug() << "end TSubLattice<T>::resetParticleRotation()\n";
}

/*!
  Send data of owned particles to the master
*/
template <class T>
void TSubLattice<T>::sendDataToMaster()
{
  console.Debug() << "TSubLattice<T>::sendDataToMaster()\n";
  vector<Vec3> positions;
  vector<double> radii;
  vector<int> ids;

  m_ppa->forAllParticlesGet(positions,(Vec3 (T::*)() const)(&T::getPos));
  m_ppa->forAllParticlesGet(radii,(double (T::*)() const)(&T::getRad));
  m_ppa->forAllParticlesGet(ids,(int (T::*)() const)(&T::getID));

  m_tml_comm.send_gather(positions,0);
  m_tml_comm.send_gather(radii,0);
  m_tml_comm.send_gather(ids,0);

  console.Debug() << "end TSubLattice<T>::sendDataToMaster()\n";
}

/*!
  Send number of owned particles to the master
*/
template <class T>
void TSubLattice<T>::countParticles()
{
  console.Debug()<<"TSubLattice<T>::countParticles()\n";
  CMPIVarSGBufferLeaf buffer(m_comm,0);
  //-- pack particles into buffer
  buffer.append(m_ppa->size());
  // send
  buffer.send();
}

/*!
  Print structural information
*/
template <class T>
void TSubLattice<T>::printStruct()
{
  cout<< "My Rank : " << m_rank << "\n" ;
  if(m_ppa!=NULL)
  {
          cout << *m_ppa << endl;
  }
}

template <class T>
void TSubLattice<T>::printData()
{
  cout << "Data: my rank : " << m_rank << "particles : \n" ;
  m_ppa->forAllParticles((void (T::*)())(&T::print));
}

template <class T>
void TSubLattice<T>::printTimes()
{
  console.Debug() << "time spent calculating force : " << forcetime << " sec\n";
  console.Debug() << "time spent communicating     : " << commtime << " sec\n";
  console.Debug() << "time spent packing           : " << packtime << " sec\n";
  console.Debug() << "time spent unpacking         : " << unpacktime << " sec\n";
}

//-------------- FIELD FUNCTIONS ----------------



/*!
  add scalar per-particle field to saver list
*/
template <class T>
void TSubLattice<T>::addScalarParticleField()
{
  // cout << "TSubLattice<T>::addScalarParticleField\n";
  string fieldname;
  int id,is_tagged;

  m_tml_comm.recv_broadcast_cont(fieldname,0);
  //cout << "recvd. fieldname: " << fieldname << "\n";
  m_tml_comm.recv_broadcast(id,0);
  //cout << "recvd. id: " << id << "\n";
  m_tml_comm.recv_broadcast(is_tagged,0);
  //cout << "recvd. is_tagged: " << is_tagged << "\n";

  typename T::ScalarFieldFunction rdf=T::getScalarFieldFunction(fieldname);
  ScalarParticleFieldSlave<T> *new_spfs;
  if(is_tagged==0)
  {
    new_spfs=new ScalarParticleFieldSlave<T>(&m_tml_comm,m_ppa,rdf);
  }
  else
  {
    int tag,mask;
    m_tml_comm.recv_broadcast(tag,0);
    console.XDebug() << "recvd. tag: " << tag << "\n";
    m_tml_comm.recv_broadcast(mask,0);
    console.XDebug() << "recvd. mask: " << mask << "\n";
    new_spfs=new ScalarParticleFieldSlaveTagged<T>(&m_tml_comm,m_ppa,rdf,tag,mask);
  }
  m_field_slaves.insert(make_pair(id,new_spfs));
}

/*!
  add vector per-particle field to saver list
*/
template <class T>
void TSubLattice<T>::addVectorParticleField()
{
  console.XDebug() << "TSubLattice<T>::addVectorParticleField\n";
  string fieldname;
  int id,is_tagged;

  m_tml_comm.recv_broadcast_cont(fieldname,0);
  console.XDebug() << "recvd. fieldname: " << fieldname << "\n";
  m_tml_comm.recv_broadcast(id,0);
  console.XDebug() << "recvd. id: " << id << "\n";
  m_tml_comm.recv_broadcast(is_tagged,0);
  console.XDebug() << "recvd. is_tagged: " << is_tagged << "\n";

  typename T::VectorFieldFunction rdf=T::getVectorFieldFunction(fieldname);
  VectorParticleFieldSlave<T> *new_vpfs;
  if(is_tagged==0)
  {
    new_vpfs=new VectorParticleFieldSlave<T>(&m_tml_comm,m_ppa,rdf);
  }
  else
  {
    int tag,mask;
    m_tml_comm.recv_broadcast(tag,0);
    console.XDebug() << "recvd. tag: " << tag << "\n";
    m_tml_comm.recv_broadcast(mask,0);
    console.XDebug() << "recvd. mask: " << mask << "\n";
    new_vpfs=new VectorParticleFieldSlaveTagged<T>(&m_tml_comm,m_ppa,rdf,tag,mask);
  }
  m_field_slaves.insert(make_pair(id,new_vpfs));

  console.Debug() << "end TSubLattice<T>::addVectorParticleField\n";
}


/*!
  add scalar per-interaction field to saver list
*/
template <class T>
void TSubLattice<T>::addScalarInteractionField()
{
  console.XDebug() << "TSubLattice<T>::addScalarInteractionField\n";
  string fieldname;
  string igname;
  string igtype;
  int id,is_tagged,tag,mask,is_checked;

  m_tml_comm.recv_broadcast_cont(fieldname,0);
  console.XDebug() << "recvd. fieldname: " << fieldname << "\n";
  m_tml_comm.recv_broadcast(id,0);
  console.XDebug() << "recvd. id: " << id << "\n";
  m_tml_comm.recv_broadcast_cont(igname,0);
  console.XDebug() << "recvd. interaction group name: " << igname << "\n";
  m_tml_comm.recv_broadcast_cont(igtype,0);
  console.XDebug() << "recvd. interaction group name: " << igtype << "\n";
  m_tml_comm.recv_broadcast(is_tagged,0);
  console.XDebug() << "recvd. is_tagged: " << is_tagged << "\n";

  // get interaction group
  map<string,AParallelInteractionStorage*>::iterator it=m_dpis.find(igname);
  if(is_tagged==1)
  {
    m_tml_comm.recv_broadcast(tag,0);
    m_tml_comm.recv_broadcast(mask,0);
  }
  m_tml_comm.recv_broadcast(is_checked,0);
  console.XDebug() << "recvd. is_checked: " << is_checked << "\n";

  if(it!=m_dpis.end())
  {
    console.XDebug() << "found interaction group in dynamic\n";
    AFieldSlave* new_sifs;
    new_sifs=(it->second)->generateNewScalarFieldSlave(&m_tml_comm,fieldname,is_checked,is_tagged,tag,mask);
    m_field_slaves.insert(make_pair(id,new_sifs));
  }
  else
  {
    it=m_bpis.find(igname);
    if(it!=m_bpis.end()){
       console.XDebug() << "found interaction group in bonded\n";
      AFieldSlave* new_sifs;
      new_sifs=(it->second)->generateNewScalarFieldSlave(&m_tml_comm,fieldname,is_checked,is_tagged,tag,mask);
      m_field_slaves.insert(make_pair(id,new_sifs));
    }
    else // not in dynamic or bonded -> try damping
      {
	//typename map<string,CDampingGroup<T>*>::iterator it2=m_damping.find(igname);
	it=m_damping.find(igname);
	if(it!=m_damping.end()) // found it in damping
	  {
	    AFieldSlave* new_sifs;
	    new_sifs=(it->second)->generateNewScalarFieldSlave(&m_tml_comm,fieldname,is_checked,is_tagged,tag,mask);
	    m_field_slaves.insert(make_pair(id,new_sifs));
	  }
	else // still not found -> unknown name -> error
	  {
            cerr << "ERROR : unknown interaction group name in addScalarInteractionField " << endl;
	  }
      }
  }

   console.XDebug() << "end TSubLattice<T>::addScalarInteractionField\n";
}

/*!
  add per-interaction vector field to saver list
*/
template <class T>
void TSubLattice<T>::addVectorInteractionField()
{
  console.Debug() << "TSubLattice<T>::addVectorInteractionField\n";
  string fieldname;
  string igname;
  string igtype;
  int id,is_tagged,tag,mask,is_checked;

  m_tml_comm.recv_broadcast_cont(fieldname,0);
  console.XDebug() << "recvd. fieldname: " << fieldname << "\n";
  m_tml_comm.recv_broadcast(id,0);
  console.XDebug() << "recvd. id: " << id << "\n";
  m_tml_comm.recv_broadcast_cont(igname,0);
  console.XDebug() << "recvd. interaction group name: " << igname << "\n";
  m_tml_comm.recv_broadcast_cont(igtype,0);
  console.XDebug() << "recvd. interaction group type: " << igtype << "\n";
  m_tml_comm.recv_broadcast(is_tagged,0);
  console.XDebug() << "recvd. is_tagged: " << is_tagged << "\n";

  // get interaction group
  map<string,AParallelInteractionStorage*>::iterator it=m_dpis.find(igname);
  if(is_tagged==1)
  {
    m_tml_comm.recv_broadcast(tag,0);
    m_tml_comm.recv_broadcast(mask,0);
  }
  m_tml_comm.recv_broadcast(is_checked,0);
  console.XDebug() << "recvd. is_checked: " << is_checked << "\n";

  if(it!=m_dpis.end())
  {
    console.XDebug() << "found interaction group in dynamic\n";
    AFieldSlave* new_sifs;
    new_sifs=(it->second)->generateNewVectorFieldSlave(&m_tml_comm,fieldname,is_checked,is_tagged,tag,mask);
    if(new_sifs!=NULL){
      m_field_slaves.insert(make_pair(id,new_sifs));
    } else {
      console.Error()<<"ERROR: could not generate Field Slave for field " << fieldname << "\n";
    }
  }
  else
  {
    it=m_bpis.find(igname);
    if(it!=m_bpis.end()){
      console.XDebug() << "found interaction group in bonded\n";
      AFieldSlave* new_sifs;
      new_sifs=(it->second)->generateNewVectorFieldSlave(&m_tml_comm,fieldname,is_checked,is_tagged,tag,mask);
      m_field_slaves.insert(make_pair(id,new_sifs));
    }
    else // not in dynamic or bonded -> try damping
      {
	//typename map<string,CDampingGroup<T>*>::iterator it2=m_damping.find(igname);
	it=m_damping.find(igname);
	if(it!=m_damping.end()) // found it in damping
	  {
	    AFieldSlave* new_sifs;
	    new_sifs=(it->second)->generateNewVectorFieldSlave(&m_tml_comm,fieldname,is_checked,is_tagged,tag,mask);
	    m_field_slaves.insert(make_pair(id,new_sifs));
	  }
	else // still not found -> unknown name -> error
	  {
            cerr << "ERROR : unknown interaction group name in addScalarInteractionField " << endl;
	  }
      }
  }

  console.Debug() << "end TSubLattice<T>::addVectorInteractionField\n";
}

/*!
  add scalar per-interaction "history" field to saver list
*/
template <class T>
void TSubLattice<T>::addScalarHistoryInteractionField()
{
  console.XDebug() << "TSubLattice<T>::addScalarHistoryInteractionField\n";
  string fieldname;
  string igname;
  string igtype;
  int id,is_tagged,tag,mask,is_checked;

  m_tml_comm.recv_broadcast_cont(fieldname,0);
  console.XDebug() << "recvd. fieldname: " << fieldname << "\n";
  m_tml_comm.recv_broadcast(id,0);
  console.XDebug() << "recvd. id: " << id << "\n";
  m_tml_comm.recv_broadcast_cont(igname,0);
  console.XDebug() << "recvd. interaction group name: " << igname << "\n";
  m_tml_comm.recv_broadcast_cont(igtype,0);
  console.XDebug() << "recvd. interaction group name: " << igtype << "\n";
  m_tml_comm.recv_broadcast(is_tagged,0);
  console.XDebug() << "recvd. is_tagged: " << is_tagged << "\n";

  if(is_tagged==1)
  {
    m_tml_comm.recv_broadcast(tag,0);
    m_tml_comm.recv_broadcast(mask,0);
  }
  m_tml_comm.recv_broadcast(is_checked,0);
  console.XDebug() << "recvd. is_checked: " << is_checked << "\n";

  // get interaction group
  map<string,AParallelInteractionStorage*>::iterator it=m_bpis.find(igname);
  if(it!=m_bpis.end()){
      console.XDebug() << "found interaction group in bonded\n";
      AFieldSlave* new_sifs;
      new_sifs=(it->second)->generateNewScalarHistoryFieldSlave(&m_tml_comm,fieldname,is_tagged,tag,mask);
      if(new_sifs!=NULL){
	m_field_slaves.insert(make_pair(id,new_sifs));
      }
  } else { // not in dynamic or bonded -> throw error
      cerr << "ERROR : unknown interaction group name in addScalarHistoryInteractionField " << endl;
  }

   console.XDebug() << "end TSubLattice<T>::addScalarHistoryInteractionField\n";
}


/*!
  Add a per-triangle (in case of a TriMesh) of per-edge (in case of a Mesh2D) vector
  field saver. Data received from master.
*/
template <class T>
void TSubLattice<T>::addVectorTriangleField()
{
  console.Debug() << "TSubLattice<T>::addVectorTriangleField()\n";
  string fieldname;
  string meshname;
  int id;

  // receive params
  m_tml_comm.recv_broadcast_cont(fieldname,0);
  console.XDebug() << "recvd. fieldname: " << fieldname << "\n";
  m_tml_comm.recv_broadcast_cont(meshname,0);
  console.XDebug() << "recvd. meshname: " << meshname << "\n";
  m_tml_comm.recv_broadcast(id,0);
  console.XDebug() << "recvd. id: " << id << "\n";

  map<string,TriMesh*>::iterator tm=m_mesh.find(meshname);
  // if meshname is in trimesh map
  if (tm!=m_mesh.end()){
    // get reader function
    Triangle::VectorFieldFunction rdf=Triangle::getVectorFieldFunction(fieldname);
    // check it
    if(rdf==NULL){
      console.Critical() << "NULL rdf for field " << fieldname << "in mesh " << meshname << "\n";
    }
    VectorTriangleFieldSlave* new_vfs=new VectorTriangleFieldSlave(&m_tml_comm,tm->second,rdf);
    m_field_slaves.insert(make_pair(id,new_vfs));
  } else {
    map<string,Mesh2D*>::iterator m2d=m_mesh2d.find(meshname);
    if(m2d!=m_mesh2d.end()){
      Edge2D::VectorFieldFunction rdf=Edge2D::getVectorFieldFunction(fieldname);
      // check it
      if(rdf==NULL){
		console.Critical() << "NULL rdf for field " << fieldname << "in mesh " << meshname << "\n";
      }
      VectorEdge2DFieldSlave* new_efs= new VectorEdge2DFieldSlave(&m_tml_comm,m2d->second,rdf);
      m_field_slaves.insert(make_pair(id,new_efs));
    } else {
		console.Critical() << "trying to setup VectorTriangleField " << fieldname << " for non-existing mesh " << meshname << "\n";
	}
  }
  console.Debug() << "end TSubLattice<T>::addVectorTriangleField()\n";
}

/*!
  Add a per-triangle scalar field saver. Data received from master.
*/
template <class T>
void TSubLattice<T>::addScalarTriangleField()
{
  console.Debug() << "TSubLattice<T>::addScalarTriangleField()\n";
  string fieldname;
  string meshname;
  int id;

  // receive params
  m_tml_comm.recv_broadcast_cont(fieldname,0);
  console.XDebug() << "recvd. fieldname: " << fieldname << "\n";
  m_tml_comm.recv_broadcast_cont(meshname,0);
  console.XDebug() << "recvd. meshname: " << meshname << "\n";
  m_tml_comm.recv_broadcast(id,0);
  console.XDebug() << "recvd. id: " << id << "\n";

  // get reader function
  Triangle::ScalarFieldFunction rdf=Triangle::getScalarFieldFunction(fieldname);
  // check it
  if(rdf==NULL){
    console.Critical() << "NULL rdf for field " << fieldname << "in mesh " << meshname << "\n";
  }
  
  map<string,TriMesh*>::iterator tm=m_mesh.find(meshname);
  // if meshname is in trimesh map
  if (tm!=m_mesh.end()){
    // get reader function
	ScalarTriangleFieldSlave* new_vtfs=new ScalarTriangleFieldSlave(&m_tml_comm,tm->second,rdf);
	m_field_slaves.insert(make_pair(id,new_vtfs));
  } else {
	console.Critical() << "trying to setup ScalarTriangleField " << fieldname << " for non-existing mesh " << meshname << "\n";
  }
  console.Debug() << "end TSubLattice<T>::addScalarTriangleField()\n";
}

/*!
  Add vector wall field. Data received from master.
*/
template <class T>
void TSubLattice<T>::addVectorWallField()
{
  console.XDebug() << "begin TSubLattice<T>::addVectorWallField()\n";
  string fieldname;
  string tmp_wallname;
  vector<string> wallnames;
  int nwall;
  int id;

  // receive params
  m_tml_comm.recv_broadcast_cont(fieldname,0);
  console.XDebug() << "recvd. fieldname: " << fieldname << "\n";
  m_tml_comm.recv_broadcast(nwall,0);
  console.XDebug() << "recvd. nwall: " << nwall << "\n";
  for(int i=0;i<nwall;i++){
    m_tml_comm.recv_broadcast_cont(tmp_wallname,0);
    wallnames.push_back(tmp_wallname);
    console.XDebug() << "recvd. wallname: " << tmp_wallname << "\n";
    tmp_wallname.clear();
  }
  m_tml_comm.recv_broadcast(id,0);
  console.XDebug() << "recvd. id: " << id << "\n";

  // check validity of 1st wall name
  map<string,CWall*>::iterator cwalliter=m_walls.find(*(wallnames.begin()));
  if(cwalliter==m_walls.end()){ // 1st wall name invalid -> exit
    std::stringstream msg;
    msg
      << "ERROR in addVectorWallField: wallname '"
      << *(wallnames.begin()) << " 'invalid. Valid wall names: ";
    for (map<string,CWall*>::const_iterator it = m_walls.begin(); it != m_walls.end(); it++)
    {
      msg << "'" << it->first << "' ";
    }
    console.Error() << msg.str() << "\n";
    throw std::runtime_error(msg.str());
  } else { // first wall valid -> try to get slave
    // get summation flag from wall
    int sumflag=(cwalliter->second)->getFieldSummationFlag(fieldname);
    // if process 1, send summation flag back to master
    if(m_tml_comm.rank()==1){
      m_tml_comm.send(sumflag,0);
    }
    m_tml_comm.barrier();

    //field slave
    AWallFieldSlave* new_fs=(cwalliter->second)->generateVectorFieldSlave(&m_tml_comm,fieldname);

    // try to insert other walls
    vector<string>::iterator niter=wallnames.begin();
    if(niter!=wallnames.end()) niter++ ; // jump past 1st wall - already got it
    while(niter!=wallnames.end()){
      string wname=*niter;
      map<string,CWall*>::iterator iter=m_walls.find(wname);
      if(iter==m_walls.end()){ // wall name invalid -> exit
        std::stringstream msg;
        msg
          << "ERROR in addVectorWallField: wallname '"
          << wname << " 'invalid";
        for (map<string,CWall*>::const_iterator it = m_walls.begin(); it != m_walls.end(); it++)
        {
          msg << "'" << it->first << "' ";
        }

        console.Error() << msg.str() << "\n";
        throw std::runtime_error(msg.str());
      } else {
      	new_fs->addWall(iter->second);
      }
      niter++;
    }
    if(new_fs!=NULL){
      m_field_slaves.insert(make_pair(id,new_fs));
    } else {
      console.Error() << "ERROR in addVectorWallField: got NULL Slave\n";
    }
  }

  console.XDebug() << "end TSubLattice<T>::addVectorWallField()\n";
}

/*!
  send field data to master
*/
template <class T>
void TSubLattice<T>::sendFieldData()
{
  console.Debug() << "TSubLattice<T>::sendFieldData()\n";
  // receive id's of field to send
  int id;
  m_tml_comm.recv_broadcast(id,0);
  console.Debug()  << "received field id " << id << " for data collection\n" ;
  if(m_field_slaves[id] != NULL)
  {
    m_field_slaves[id]->sendData();
  }
  else
  { // can not happen
    cerr << "NULL pointer in m_field_slaves!" << endl;
  }
  // call the sending function of the field
  console.Debug() << "end TSubLattice<T>::sendFieldData()\n";
}


// ---- Checkpointing ----------
/*!
  save snapshot data, i.e. for viz/postprocessing
*/
template <class T>
void TSubLattice<T>::saveSnapShotData(std::ostream &oStream)
{
  // get precision of output stream and set it to 9 significant digits
  std::streamsize oldprec=oStream.precision(9);

  //
  // Save particle check-point data
  //
  ParticleArray &particleArray = dynamic_cast<ParticleArray &>(*m_ppa);
  typename ParticleArray::ParticleIterator
    particleIt(particleArray.getInnerParticleIterator());

  const std::string delim = "\n";

  oStream << particleIt.getNumRemaining() << delim;
  while (particleIt.hasNext()) {
    particleIt.next().saveSnapShotData(oStream);
    oStream << delim;
  }

  //
  // Save Bonded interaction check-point data.
  //
  typedef std::map<string,AParallelInteractionStorage*> NameBondedInteractionsMap;
  typename NameBondedInteractionsMap::iterator it;
  oStream << m_bpis.size() << delim;
  for (it = m_bpis.begin(); it != m_bpis.end(); it++) {
    it->second->saveSnapShotData(oStream);
    oStream << delim;
  }

  // dump trimeshdata (if exists)
  oStream << "TMIG " << m_mesh.size() << delim;
  for(typename map<string,TriMesh*>::iterator tm_iter=m_mesh.begin();
      tm_iter!=m_mesh.end();
      tm_iter++){
    oStream << tm_iter->first << delim;
    tm_iter->second->writeCheckPoint(oStream,delim);
  }

  // restore output stream to old precision
  oStream.precision(oldprec);
}

/*!
  save checkpoint data, i.e. for restarting
*/
template <class T>
void TSubLattice<T>::saveCheckPointData(std::ostream &oStream)
{
  const std::string delim = "\n";
  //
  // Save particle check-point data
  //
  m_ppa->saveCheckPointData(oStream);

  //
  // Save Bonded interaction check-point data.
  //
  typedef std::map<string,AParallelInteractionStorage*> NameBondedInteractionsMap;
  typename NameBondedInteractionsMap::iterator it;
  oStream << m_bpis.size() << delim;
  for (it = m_bpis.begin(); it != m_bpis.end(); it++) {
    it->second->saveCheckPointData(oStream);
    oStream << delim;
  }

  //
  // Save Non-bonded interaction check-point data
  //
  int count_save=0;
  for(std::map<string,AParallelInteractionStorage*>::iterator iter=m_dpis.begin();
      iter!=m_dpis.end();
      iter++){
    if(iter->second->willSave()) count_save++;
  }
  oStream << count_save << delim;
  for(std::map<string,AParallelInteractionStorage*>::iterator iter=m_dpis.begin();
      iter!=m_dpis.end();
      iter++){
    if(iter->second->willSave()) {
      iter->second->saveCheckPointData(oStream);
      oStream << delim;
    }
  }

  // Save walls (name, pos, normal)
  oStream << "Walls " << m_walls.size() << delim;
  for(map<string,CWall*>::iterator w_iter=m_walls.begin();
      w_iter!=m_walls.end();
      w_iter++){
    oStream << w_iter->first << delim;
    w_iter->second->writeCheckPoint(oStream,delim);
  }

  // Save sphere bodies (name, pos, radius)
  oStream << "Spheres " << m_spheres.size() << delim;
  for(map<string,CSphereBody*>::iterator w_iter=m_spheres.begin();
      w_iter!=m_spheres.end();
      w_iter++){
    oStream << w_iter->first << delim;
    w_iter->second->writeCheckPoint(oStream,delim);
  }

  // dump trimeshdata (if exists)
  oStream << "TriMesh " << m_mesh.size() << delim;
  for(typename map<string,TriMesh*>::iterator tm_iter=m_mesh.begin();
      tm_iter!=m_mesh.end();
      tm_iter++){
    oStream << tm_iter->first << delim;
    tm_iter->second->writeCheckPoint(oStream,delim);
  }
  // dump 2D mesh data (if exists)
  oStream << "Mesh2D " << m_mesh2d.size() << delim;
  for(typename map<string,Mesh2D*>::iterator tm_iter=m_mesh2d.begin();
      tm_iter!=m_mesh2d.end();
      tm_iter++){
    oStream << tm_iter->first << delim;
    tm_iter->second->writeCheckPoint(oStream,delim);
  }
}

template <class T>
void TSubLattice<T>::loadCheckPointData(std::istream &iStream)
{
  // get particles
  m_ppa->loadCheckPointData(iStream);

  // rebuild neighbor table
  CMPIBarrier barrier(getWorkerComm());
  m_ppa->rebuild();
  barrier.wait("PPA rebuild");

  //-- get bonds --
  // get nr. of bonded interaction group in the checkpoint file
  unsigned int nr_bonded_ig;
  iStream >> nr_bonded_ig;

  // compare with existing bonded particle interaction storage (bpis)
  // barf if not equal
  if(nr_bonded_ig!=m_bpis.size()){
    std::cerr << "number of bonded interaction groups differ between snapshot and script!" << std::endl;
  } else { // numbers fit -> read data
    for (map<string,AParallelInteractionStorage*>::iterator it = m_bpis.begin();
	 it != m_bpis.end();
	 it++) { // for all interaction groups
      it->second->loadCheckPointData(iStream);
    }
  }
  //-- get nonbonded interactions --
  // get nr. of non-bonded interaction group in the checkpoint file
  unsigned int nr_nonbonded_ig;
  iStream >> nr_nonbonded_ig;

  // compare with existing non-bonded (dynamic) particle interaction storage (dpis)
  // barf if not equal
  if(nr_nonbonded_ig!=m_dpis.size()){
    std::cerr << "number of dynamic interaction groups differ between snapshot and script!" << std::endl;
  } else { // numbers fit -> read data
    for (map<string,AParallelInteractionStorage*>::iterator it = m_dpis.begin();
	 it != m_dpis.end();
	 it++) { // for all interaction groups
      it->second->loadCheckPointData(iStream);
    }
  }
  //--- load walls ---
  string token;
  iStream >> token;
  if(token!="Walls") { // found wrong token -> barf
    std::cerr << "expected Walls , got " << token << std::endl;
  }
  // nr. of walls
  int nwalls;
  iStream >> nwalls;
  // read wall names & data
  string wname;
  for(int i=0;i<nwalls;i++){
    CWall* new_wall = new CWall();
    iStream >> wname;
    new_wall->loadCheckPoint(iStream);
    m_walls[wname]=new_wall;
  }

  iStream >> token;
  if(token!="Spheres") { // found wrong token -> barf
    std::cerr << "expected Spheres , got " << token << std::endl;
  }

  // nr. of sphere bodies
  int nspheres;
  iStream >> nspheres;
  // read sphere body names & data
  string sname;
  for(int i=0;i<nspheres;i++){
    CSphereBody* new_sphere = new CSphereBody();
    iStream >> sname;
    new_sphere->loadCheckPoint(iStream);
    m_spheres[sname]=new_sphere;
  }
  // --- load meshes --
  int nmesh;
  string mname;
  // Trimesh (3D)
  iStream >> token;
  if(token!="TriMesh") { // found wrong token -> barf
    std::cerr << "expected TriMesh , got " << token << std::endl;
  }
  // nr. of meshes
  iStream >> nmesh;
  // read wall names & data
  for(int i=0;i<nmesh;i++){
    TriMesh* new_tm=new TriMesh();
    iStream >> mname;
    new_tm->loadCheckPoint(iStream);
    m_mesh.insert(make_pair(mname,new_tm));
  }
  // Mesh2D
  iStream >> token;
  if(token!="Mesh2D") { // found wrong token -> barf
    std::cerr << "expected Mesh2D , got " << token << std::endl;
  }
  // nr. of meshes
  iStream >> nmesh;
  // read wall names & data
  for(int i=0;i<nmesh;i++){
    Mesh2D* new_m2d=new Mesh2D();
    iStream >> mname;
    new_m2d->loadCheckPoint(iStream);
    m_mesh2d.insert(make_pair(mname,new_m2d));
  }
}

// -- mesh data exchange functions --

/*!
  send back node ids of given mesh, get mesh name from master
*/
template <class T>
void TSubLattice<T>::getMeshNodeRef()
{
  console.XDebug() << "TSubLattice<T>::getMeshNodeRef()\n";
  vector<int> ref_vec;

  // MPI buffer
  CVarMPIBuffer param_buffer(m_comm);
  // receive mesh name
  param_buffer.receiveBroadcast(0);
  string meshname=param_buffer.pop_string();
  console.XDebug() << "Mesh name: " << meshname << "\n";

  // find mesh & collect node ids into array
  map<string,TriMesh*>::iterator tm=m_mesh.find(meshname);
  if (tm!=m_mesh.end()){
    for(TriMesh::corner_iterator iter=(tm->second)->corners_begin();
	iter!=(tm->second)->corners_end();
	iter++){
      ref_vec.push_back(iter->getID());
    }
  } else {
    map<string,Mesh2D*>::iterator m2d=m_mesh2d.find(meshname);
    if(m2d!=m_mesh2d.end()){
      for(Mesh2D::corner_iterator iter=(m2d->second)->corners_begin();
	  iter!=(m2d->second)->corners_end();
	  iter++){
	ref_vec.push_back(iter->getID());
      }
    } else {
      console.Critical() << "ERROR - WRONG MESH NAME IN getMeshNodeRef() !! \n";
    }
  }
  // send back to master
  m_tml_comm.send_gather(ref_vec,0);

  console.XDebug() << "end TSubLattice<T>::getMeshNodeRef()\n";
}

/*!
  send back face (edge in 2D, triangle in 3d) ids of given mesh, get mesh name from master
*/
template <class T>
void TSubLattice<T>::getMeshFaceRef()
{
  console.XDebug() << "TSubLattice<T>::getMeshFaceRef()\n";
  vector<int> ref_vec;

  // MPI buffer
  CVarMPIBuffer param_buffer(m_comm);
  // receive mesh name
  param_buffer.receiveBroadcast(0);
  string meshname=param_buffer.pop_string();
  console.XDebug() << "Mesh name: " << meshname << "\n";

  // find mesh & collect node ids into array
  map<string,TriMesh*>::iterator tm=m_mesh.find(meshname);
  if (tm!=m_mesh.end()){
    for(TriMesh::triangle_iterator iter=(tm->second)->triangles_begin();
	iter!=(tm->second)->triangles_end();
	iter++){
      ref_vec.push_back(iter->getID());
    }
  } else {
    map<string,Mesh2D*>::iterator m2d=m_mesh2d.find(meshname);
    if(m2d!=m_mesh2d.end()){
      for(Mesh2D::edge_iterator iter=(m2d->second)->edges_begin();
	  iter!=(m2d->second)->edges_end();
	  iter++){
	ref_vec.push_back(iter->getID());
      }
    } else {
      console.Critical() << "ERROR - WRONG MESH NAME IN getMeshNodeRef() !! \n";
    }
  }
  // send back to master
  m_tml_comm.send_gather(ref_vec,0);

  console.XDebug() << "end TSubLattice<T>::getMeshNodeRef()\n";
}

/*!
  send back stress on faces of  given mesh, get mesh name from master
*/
template <class T>
void TSubLattice<T>::getMesh2DStress()
{
  console.XDebug() << "TSubLattice<T>::getMesh2DStress()\n";
  // receive mesh name
  // MPI buffer
  CVarMPIBuffer param_buffer(m_comm);
  param_buffer.receiveBroadcast(0);
  string meshname=param_buffer.pop_string();
  console.XDebug() << "Mesh name: " << meshname << "\n";

  // find mesh & collect data
  map<string,Mesh2D*>::iterator m2d=m_mesh2d.find(meshname);
  if(m2d!=m_mesh2d.end()){
    vector<pair<int,Vec3> > data_vec;
    // get data
    data_vec=(m2d->second)->forAllEdgesGetIndexed(&Edge2D::getForceDensity);

    // send data to master
    m_tml_comm.send_gather(data_vec,0);
  } else {
    console.Critical() << "ERROR - WRONG MESH NAME IN getMesh2DStress() !! \n";
  }

  console.XDebug() << "end TSubLattice<T>::getMesh2DStress()\n";
}

/*!
  send back stress on faces of  given mesh, get mesh name from master
*/
template <class T>
void TSubLattice<T>::getTriMeshForce()
{
  console.XDebug() << "TSubLattice<T>::getTriMeshStress(): enter\n";
  // receive mesh name
  // MPI buffers
  CVarMPIBuffer param_buffer(m_comm);
  param_buffer.receiveBroadcast(0);
  const std::string meshName = param_buffer.pop_string();
  console.XDebug() << "Mesh name: " << meshName << "\n";

  // find mesh & collect data
  map<string,TriMesh*>::iterator m=m_mesh.find(meshName);
  if(m != m_mesh.end()){
    vector<pair<int,Vec3> > data_vec;
    // get data
    data_vec=(m->second)->forAllTrianglesGetIndexed(&Triangle::getForce);

    // send data to master
    m_tml_comm.send_gather(data_vec,0);
  } else {
    std::stringstream msg;
    msg << "Invalid mesh name: " << meshName << ". No such triangular mesh.";
    throw std::runtime_error(msg.str().c_str());
  }

  console.XDebug() << "TSubLattice<T>::getTriMeshStress(): exit\n";
} 

/*!
*/
template <class T>
void TSubLattice<T>::setInteractionParameter()
{
	console.XDebug() << "TSubLattice<T>::setInteractionParameter(): enter\n";
	// receive and unpack parameters
	CVarMPIBuffer buffer(m_comm);
	buffer.receiveBroadcast(0);
	const std::string igname(buffer.pop_string());
	const std::string pname(buffer.pop_string());
	double val = buffer.pop_double();
	
	// debug output
	console.XDebug() << "received: " << igname << " , " << pname << " , " << val << "\n";
	
	// get interaction group
	map<string,AParallelInteractionStorage*>::iterator it=m_dpis.find(igname);
	if(it!=m_dpis.end()){ // found in dynamic interaction groups (elastic, friction etc...)
		console.XDebug() << "found interaction group in dynamic\n";
		(it->second)->forAllInteractionsSet(pname,val);
	} else {
		console.Error() << "TSubLattice<T>::setInteractionParameter() - interaction group " << igname << " not found\n";
	}
	
	console.XDebug() << "TSubLattice<T>::setInteractionParameter(): exit\n";
}

