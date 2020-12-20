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

#include "Parallel/SubLatticeControler.h"

// -- project includes --
#include "Parallel/mpi_tag_defs.h"
#include "Parallel/sublattice_cmd.h"
#include "Parallel/mpicmdbuf.h"
#include "Parallel/mpibarrier.h"
#include "Foundation/console.h"
#include "Model/Particle.h"
#include "Model/RotParticle.h"
#include "Model/RotParticleVi.h"
#include "Model/RotThermParticle.h"
#include "Parallel/CheckPointer.h"
#include "Parallel/RotSubLattice.h"

#include "Parallel/MpiWrap.h"

#include <stdexcept>

CSubLatticeControler::CSubLatticeControler()
  : m_lattice(NULL),
    m_pCheckPointer(NULL),
    m_pSnapShooter(NULL),
    m_timingFileName(),
    m_timersPtr(new MpiWTimers())
{
  //   console.Debug() << "CSubLatticeControler::CSubLatticeControler()" << "\n";

  m_pCheckPointer = new CheckPointer(*this);
  m_pSnapShooter = new CheckPointer(*this);
  m_timersPtr->zeroise("RebuildParticleArray");
  m_timersPtr->zeroise("RebuildInteractions");
  m_timersPtr->zeroise("NeighbourSearch");
  m_timersPtr->zeroise("CheckNeighbours");
  m_timersPtr->zeroise("CheckPoint");
  m_timersPtr->zeroise("BoundaryDataExchange");
  m_timersPtr->zeroise("ForceCalculation");
  m_timersPtr->zeroise("UpdateInteractions");
  m_timersPtr->zeroise("UpdateBondedInteractions");
  m_timersPtr->zeroise("UpdateDynamicInteractions");
}

CSubLatticeControler::~CSubLatticeControler()
{
  console.Debug()
    << "CSubLatticeControler::~CSubLatticeControler(): slave "
    << m_global_rank << " finalizing...\n";
  console.Debug()
    << "CSubLatticeControler::~CSubLatticeControler():"
    << "deleting checkpointer...\n";
  delete m_pCheckPointer;
  delete m_pSnapShooter;
  console.Debug()
    << "CSubLatticeControler::~CSubLatticeControler():"
    << "deleting sublattice...\n";
  delete m_lattice;
  // cleanup MPI Groups & Communicators
  MPI_Group_free(&m_global_group);
  MPI_Group_free(&m_local_group);

  MPI_Comm_free(&m_local_comm);

  console.Debug()
    << "CSubLatticeControler::~CSubLatticeControler(): exit.\n";
}


/*!
  Initialize MPI communicators
*/
void CSubLatticeControler::initMPI()
{
  std::cout << "CSubLatticeControler::initMPI()\n";

  // global communicator = MPI_COMM_WORLD
  m_global_comm=MPI_COMM_WORLD;
  m_tml_global_comm.setComm(m_global_comm);

  // -- local communicator = global comm - proc(0)
  // get global MPI_Group
  MPI_Comm_group(MPI_COMM_WORLD,&m_global_group);
  // subtract id 0 from global group
  int id0=0;
  MPI_Group_excl(m_global_group,1,&id0,&m_local_group);
  //  create communicator
  MPI_Comm_create(MPI_COMM_WORLD,m_local_group,&m_local_comm);
  m_tml_local_comm.setComm(m_local_comm);

  m_tml_global_comm.barrier("CSubLatticeControler::initMPI");

  // get size of the global communicator and my rank in it
  m_global_size=m_tml_global_comm.size();
  m_global_rank=m_tml_global_comm.rank();
  // get size of the local (i.e. workers only) communicator and my rank in it
  m_local_size=m_tml_local_comm.size();
  m_local_rank=m_tml_local_comm.rank();
  // tell checkpointers about the mpi communicators
  m_pCheckPointer->setMpiComm(m_global_comm);
  m_pSnapShooter->setMpiComm(m_global_comm);

  std::cout << "slave started at local/global rank " << m_local_rank << " / " << m_global_rank << "\n";
}

/*!
  make a new Lattice of the correct type

  \todo make it throw an exeption if it fails
*/
void CSubLatticeControler::makeLattice()
{
  console.Debug() << "CSubLatticeControler::makeLattice()\n" ;

  CVarMPIBuffer param_buffer(m_global_comm);

  // get lattice parameters
  param_buffer.receiveBroadcast(0);

  console.Debug() << "got params " << "\n";

  CLatticeParam clp = CLatticeParam::extractLatticeParam(&param_buffer);

  // make the lattice, depending of the type
  if (clp.particle_type()=="Basic") {
    console.Debug() << "making lattice with basic particles\n";
    m_lattice=new TSubLattice<CParticle>(clp,m_global_rank,m_global_comm,m_local_comm);
  } else if (clp.particle_type()=="Rot") {
    console.Debug() << "making lattice with simple rotational particles\n";
    m_lattice=new TRotSubLattice<CRotParticle>(clp,m_global_rank,m_global_comm,m_local_comm);
  } else if (clp.particle_type()=="RotVi") {
    console.Debug() << "making lattice with simple Verlet-integration rotational particles\n";
    m_lattice = new TRotSubLattice<CRotParticleVi>(clp,m_global_rank,m_global_comm,m_local_comm);
  } else if (clp.particle_type()=="RotTherm") {
    console.Debug() << "making lattice with thermal Verlet-integration rotational particles\n";
    m_lattice = new TRotSubLattice<CRotThermParticle>(clp,m_global_rank,m_global_comm,m_local_comm);
  } else {
    std::stringstream msg;
    msg
      << "ERROR! Trying to make lattice from unknown particles: "
      << clp.particle_type();
    throw std::runtime_error(msg.str());
  }
  m_lattice->setTimer(*m_timersPtr);

  console.Debug() << "end CSubLatticeControler::makeLattice()\n" ;
}

void CSubLatticeControler::saveCheckPointData(std::ostream &oStream)
{
  if (m_lattice != NULL) {
    m_lattice->saveCheckPointData(oStream);
  }
}

void CSubLatticeControler::saveSnapShotData(std::ostream &oStream)
{
  if (m_lattice != NULL) {
    m_lattice->saveSnapShotData(oStream);
  }
}

void CSubLatticeControler::loadCheckPointData(std::istream &iStream)
{
  if (m_lattice != NULL) {
    m_lattice->loadCheckPointData(iStream);
  }
}

/*!
  Initialize lattice. All boundaries are assumed to be open, i.e. not circular.
  Recieves all necessary data from Master
*/
void CSubLatticeControler::initLattice()
{
  console.XDebug() << "CSubLatticeControler::initLattice: enter\n";

  // get corners
  vector<Vec3> bbox;
  console.XDebug()
    << "CSubLatticeControler::initLattice: receiving bbox\n";

  m_tml_global_comm.recv_broadcast_cont_packed(bbox,0);
  if(bbox.size()!=2){
    console.Critical()<<"CRITICAL ERROR: bbox.size()!=2\n";
  }
  console.XDebug()
    << "CSubLatticeControler::initLattice: initialising ntable,"
    << "bbox = " << bbox[0] << ", " << bbox[1] << "\n";

  // init neighbbor table
  m_lattice->initNeighborTable(bbox[0],bbox[1]);

  // send back ppa dims and coords
  console.XDebug()
    << "CSubLatticeControler::initLattice: sending dims\n";
  vector<int> commdims=m_lattice->getCommDims();
  m_tml_global_comm.send_gather(commdims,0);

  vector<int> commcoords=m_lattice->getCommCoords();
  console.XDebug()
    << "CSubLatticeControler::initLattice: sending coords\n";
  m_tml_global_comm.send_gather(commcoords,0);

  console.XDebug() << "CSubLatticeControler::initLattice: exit\n";
}

void CSubLatticeControler::setTimeStepSize()
{
  console.Debug() << "enter CSubLatticeControler::setTimeStepSize()\n" ;
  CVarMPIBuffer buffer(m_global_comm);
  // get double parameter
  console.Debug()
    << "CSubLatticeControler::setTimeStepSize: receiving data...\n";
  buffer.receiveBroadcast(0);
  console.Debug()
    << "CSubLatticeControler::setTimeStepSize: broadcast received\n";

  console.Debug()
    << "CSubLatticeControler::setTimeStepSize: popping double...\n";
  const double dt = buffer.pop_double();
  console.Debug()
    << "CSubLatticeControler::setTimeStepSize: popped double = " << dt << "\n";

  console.Debug()
    << "CSubLatticeControler::setTimeStepSize: setting time step size\n";
  m_lattice->setTimeStepSize(dt);
  console.Debug()
    << "CSubLatticeControler::setTimeStepSize: done setting time step size.\n";

  console.Debug() << "end   CSubLatticeControler::setTimeStepSize()\n" ;
}

/*!
  Initialize lattice with at least one circular boundary condition.
  Receives all necessary data from Master
*/
void CSubLatticeControler::initLatticeCirc()
{
  console.Debug() << "CSubLatticeControler::initLatticeCirc: enter\n";

  // get corners
  vector<Vec3> bbox;
  console.Debug()
    << "CSubLatticeControler::initLatticeCirc:"
    << " receiving bounding box...\n";
  m_tml_global_comm.recv_broadcast_cont_packed(bbox,0);
  if(bbox.size()!=2){
    console.Error()<<"CRITICAL ERROR: bbox.size()!=2\n";
  }
  // get circular/non-circular boundary conditions
  vector<int> cv;
  vector<bool> circ;
  console.Debug()
    << "CSubLatticeControler::initLatticeCirc:"
    << " circular boundary prms...\n";
  m_tml_global_comm.recv_broadcast_cont(cv,0);
  for(vector<int>::iterator iter=cv.begin();
      iter!=cv.end();
      iter++){
    bool c=(*iter==1);
    console.XDebug() << "circ " << c << "\n";
    circ.push_back(c);
  }

  console.Debug()
    << "CSubLatticeControler::initLatticeCirc:"
    << " initialising neighbour table...\n";
  m_lattice->initNeighborTable(bbox[0],bbox[1],circ);

  console.Debug()
    << "CSubLatticeControler::initLatticeCirc:"
    << " sending mpi-comm dimensions...\n";
  vector<int> commdims=m_lattice->getCommDims();
  m_tml_global_comm.send_gather(commdims,0);

  console.Debug()
    << "CSubLatticeControler::initLatticeCirc:"
    << " sending mpi-comm coordinates...\n";
  vector<int> commcoords=m_lattice->getCommCoords();
  m_tml_global_comm.send_gather(commcoords,0);

  console.Debug() << "CSubLatticeControler::initLatticeCirc: exit\n";
}

void CSubLatticeControler::performTiming()
{
  CVarMPIBuffer buffer(m_global_comm);

  // get check-point parameters
  buffer.receiveBroadcast(0);
  std::string fileNamePrefix = buffer.pop_string();
  std::stringstream sStream;
  sStream << fileNamePrefix << MpiInfo(m_global_comm).rank() << ".csv";
  setTimingFileName(sStream.str());
}

void CSubLatticeControler::saveTimingData()
{
  if (getTimingFileName() != "") {
    m_timersPtr->appendData(getTimingFileName());
    m_timersPtr->zeroise();
  }
}

void CSubLatticeControler::searchNeighbors()
{
  m_lattice->searchNeighbors();
  return;
}

void CSubLatticeControler::do2dCalculations()
{
  console.Debug() << "enter CSubLatticeControler::do2dCalculations()\n" ;
  CVarMPIBuffer buffer(m_global_comm);
  CMPIBarrier barrier(m_global_comm);
  // get boolean parameter
  console.Debug()
    << "CSubLatticeControler::do2dCalculations: receiving data...\n";
  buffer.receiveBroadcast(0);
  console.Debug()
    << "CSubLatticeControler::do2dCalculations: broadcast received\n";

  console.Debug()
    << "CSubLatticeControler::do2dCalculations: popping bool...\n";
  const bool do2d = static_cast<bool>(buffer.pop_int());
  console.Debug()
    << "CSubLatticeControler::do2dCalculations: popped bool\n";

  console.Debug()
    << "CSubLatticeControler::do2dCalculations: setting bool\n";
  m_lattice->do2dCalculations(do2d);
  console.Debug()
    << "CSubLatticeControler::do2dCalculations: done setting bool.\n";

  console.Debug() << "end   CSubLatticeControler::do2dCalculations()\n" ;
}

void CSubLatticeControler::getNumParticles()
{
  int numParticles = m_lattice->getNumParticles();
  vector<int> numParticlesVector(1, numParticles);
  m_tml_global_comm.send_gather(numParticlesVector, 0);
}

void CSubLatticeControler::findParticleNearestToPoint()
{
  console.Debug()
    << "CSubLatticeControler::findParticleNearestToPoint: enter\n";

  CVarMPIBuffer buffer(m_global_comm);
  // get boolean parameter
  console.Debug()
    << "CSubLatticeControler::findParticleNearestToPoint: receiving point...\n";
  buffer.receiveBroadcast(0);
  const Vec3 pt = buffer.pop_vector();

  vector<std::pair<double,int> >
    distIdVector(1, m_lattice->findParticleNearestTo(pt));
  console.Debug()
    << "CSubLatticeControler::findParticleNearestToPoint: "
    << "sending distance data...\n";
  m_tml_global_comm.send_gather(distIdVector, 0);
  console.Debug()
    << "CSubLatticeControler::findParticleNearestToPoint: exit\n";
}

void CSubLatticeControler::getParticlePosn()
{
  console.Debug()
    << "CSubLatticeControler::getParticlePosn: enter\n";

  CVarMPIBuffer buffer(m_global_comm);
  // get boolean parameter
  console.Debug()
    << "CSubLatticeControler::getParticlePosn: receiving particle id...\n";
  buffer.receiveBroadcast(0);
  const int particleId = buffer.pop_int();

  vector<std::pair<int,Vec3> >
    idPosnVector(1, m_lattice->getParticlePosn(particleId));
  console.Debug()
    << "CSubLatticeControler::getParticlePosn: sending posn data...\n";
  m_tml_global_comm.send_gather(idPosnVector, 0);
  console.Debug()
    << "CSubLatticeControler::getParticlePosn: exit\n";
}

void CSubLatticeControler::getIdParticleData()
{
  console.Debug()
    << "CSubLatticeControler::getIdParticleData: enter\n";

  ASubLattice::IdVector particleIdVector;
  console.Debug()
    << "CSubLatticeControler::getIdParticleData: receiving particle id's\n";
  m_tml_global_comm.recv_broadcast_cont(particleIdVector, 0);
  console.Debug()
    << "CSubLatticeControler::getIdParticleData:"
    << " received " << particleIdVector.size() << " particle id's\n";

  m_lattice->getParticleData(particleIdVector);
  console.Debug()
    << "CSubLatticeControler::getIdParticleData: exit\n";
}

void CSubLatticeControler::moveSingleParticle()
{
  console.Debug() << "CSubLatticeControler::moveSingleParticle: enter\n";
  CVarMPIBuffer buffer(m_global_comm);

  console.Debug()
    << "CSubLatticeControler::moveSingleParticle:"
    " receiving id and position\n";
  buffer.receiveBroadcast(0); // get data from master
  const int particleId = buffer.pop_int();
  const Vec3 posn      = buffer.pop_vector();
  console.Debug()
    << "CSubLatticeControler::moveSingleParticle:"
    << " received id=" << particleId << " and position=" << posn << "\n";

  m_lattice->moveSingleParticleTo(particleId, posn);
  console.Debug() << "CSubLatticeControler::moveSingleParticle: exit\n";
}

/*!
  Translate mesh by given amount. Receive data from master
  and call function in SubLattice with the received parameters
*/
void CSubLatticeControler::translateMeshBy()
{
  console.Debug() << "CSubLatticeControler::translateMeshBy: enter\n";
  CVarMPIBuffer buffer(m_global_comm);

  console.Debug()
    << "CSubLatticeControler::translateMeshBy:"
    " receiving mesh-name and translation\n";
  buffer.receiveBroadcast(0); // get data from master
  const std::string meshName(buffer.pop_string());
  const Vec3 translation(buffer.pop_vector());
  console.Debug()
    << "CSubLatticeControler::translateMeshBy:"
    << " received name=\"" << meshName
    << "\" and translation=(" << translation << ")\n";

  m_lattice->translateMeshBy(meshName, translation);
  console.Debug() << "CSubLatticeControler::translateMeshBy: exit\n";
}


/*!
  Rotate mesh around axis by given amount. Receive data from master
  and call function in SubLattice with the received parameters
*/
void CSubLatticeControler::rotateMeshBy()
{
    console.Debug() << "CSubLatticeControler::rotateMeshBy: enter\n";
    CVarMPIBuffer buffer(m_global_comm);

    buffer.receiveBroadcast(0); // get data from master
    const std::string meshName(buffer.pop_string());
    const Vec3 origin(buffer.pop_vector());
    const Vec3 axis(buffer.pop_vector());
    const double angle=buffer.pop_double();
    console.Debug()
        << "CSubLatticeControler::rotateMeshBy:"
        << " received name= " << meshName
        << " origin=(" << origin << " ) axis=( " << axis << ")  angle " << angle << "\n";

    m_lattice->rotateMeshBy(meshName, origin, axis, angle);
    console.Debug() << "CSubLatticeControler::rotateMeshBy: exit\n";
}

/*!
  Set console verbosity. level recieved from master.
*/
void CSubLatticeControler::setVerbosity()
{
  CVarMPIBuffer buffer(m_global_comm);
  buffer.receiveBroadcast(0);
  const int verbose = buffer.pop_int();
  console.SetVerbose(verbose);
}

/*!
  Initialize console. Parameters recieved from master.
*/
void CSubLatticeControler::initializeConsole()
{
  CVarMPIBuffer buffer(m_global_comm);
  buffer.receiveBroadcast(0);
  const std::string filename(buffer.pop_string());
  const int bufflen = buffer.pop_int();
  console.Initialize(filename);
  console.SetBuffered(bufflen);
}

/*!
  Set console filename. filename recieved from master.
*/
void CSubLatticeControler::setConsoleFilename()
{
  CVarMPIBuffer buffer(m_global_comm);
  buffer.receiveBroadcast(0);
  const std::string filename(buffer.pop_string());
  console.SetFilename(filename);
}

/*!
  Set buffer size & mode filename. parameters recieved from master.
*/
void CSubLatticeControler::setConsoleBuffered()
{
  CVarMPIBuffer buffer(m_global_comm);
  buffer.receiveBroadcast(0);
  const int bsize = buffer.pop_int();
  console.SetBuffered(bsize);
}


/*!
  gets command messages from master and calls the functions of SubLattice
*/
void CSubLatticeControler::run()
{
  CMPILCmdBuffer cmd_buffer(m_global_comm,0);
  int command=0;
  bool is_error=false;

  #ifdef _ENABLE_DEBUG_FILE
  console.SetVerbose();
  #endif
  do{
    command=cmd_buffer.receive();
    console.Debug()
      << "CSubLatticeControler::run:"
      << "rank=" << m_global_rank
      << "->Received command=" << command << "\n";
    switch(command){
    /****fluid contents: begin****/
    case CMD_ADDFLUID : m_lattice->addFluidInteraction(); break;
    case CMD_ADD_SFF : m_lattice->addScalarFluidField(); break;
    case CMD_ADD_VFF : m_lattice->addVectorFluidField(); break;
    case CMD_ADD_SFIF : m_lattice->addScalarFluidInteractionField(); break;
    case CMD_ADD_VFIF : m_lattice->addVectorFluidInteractionField(); break;
    case CMD_SEND_COEFFI : m_lattice->sendCoeffi(); break;
    case CMD_RECV_PRESSURE : m_lattice->recvPressure(); break;
    case CMD_SOLVE: m_lattice->solveMatrix(); break;
    case CMD_UPDATE_FLUID: m_lattice->updateFluid();break;
    case CMD_EXCHANGE_CELLS: m_lattice->exchangeCells();break;
    /****fluid contents: end****/

    case CMD_PRINT :  m_lattice->printData(); break;
    case CMD_CALC :
      {
        m_timersPtr->start("ForceCalculation");
        m_lattice->oneStep();
        m_timersPtr->stop("ForceCalculation");
        break;
      }
    case CMD_XCHANGE :
      {
        m_timersPtr->start("BoundaryDataExchange");
        m_lattice->exchangePos();
        m_timersPtr->stop("BoundaryDataExchange");
        break;
      }
    case CMD_NSEARCH : searchNeighbors(); break;
    case CMD_CHECKNEIGHBORS :
      {
        m_timersPtr->start("CheckNeighbours");
        m_lattice->checkNeighbors();
        m_timersPtr->stop("CheckNeighbours");
        break;
      }
    case CMD_UPDATE :
      {
        m_timersPtr->start("UpdateInteractions");
        m_lattice->updateInteractions();
        m_timersPtr->stop("UpdateInteractions");
        break;
      }
    case CMD_CHGRADIUS : m_lattice->changeRadiusBy(); break;
    case CMD_PMOVE : m_lattice->moveParticleTo(); break;
    case CMD_PMOVETAGGEDBY : m_lattice->moveTaggedParticlesBy(); break;
    case CMD_PSETND : m_lattice->setParticleNonDynamic(); break;
    case CMD_PSETNR : m_lattice->setParticleNonRot(); break;
    case CMD_PVEL : m_lattice->setParticleVelocity(); break;
    case CMD_PANGVEL : m_lattice->setParticleAngularVelocity(); break;
    case CMD_PDENS : m_lattice->setParticleDensity(); break;
    case CMD_RPROT : m_lattice->resetParticleRotation(); break;
    case CMD_PTVEL : m_lattice->setTaggedParticleVel(); break;
    case CMD_WMOVE : m_lattice->moveWallBy(); break;
    case CMD_SPHEREBODYMOVE : m_lattice->moveSphereBodyBy(); break;
    case CMD_WNORM : m_lattice->setWallNormal(); break;
    case CMD_WFORCE : m_lattice->applyForceToWall(); break;
    case CMD_COUNT : m_lattice->countParticles(); break;
    case CMD_ADDBONDEDIG : m_lattice->addBondedIG(); break;
    case CMD_ADDSHORTBONDEDIG : m_lattice->addShortBondedIG(); break;
    case CMD_ADDCAPPEDBONDEDIG : m_lattice->addCappedBondedIG(); break;
    case CMD_ADDROTBONDEDIG : m_lattice->addRotBondedIG(); break;
    case CMD_ADDBRITTLEBEAMSCIG : m_lattice->addBrittleBeamSCIG(); break;
    case CMD_ADDBRITTLEBEAMDZCIG : m_lattice->addBrittleBeamDZCIG(); break;
    case CMD_ADDROTTHERMBONDEDIG : m_lattice->addRotThermBondedIG(); break;
    case CMD_ADDPIG : m_lattice->addPairIG(); break;
    case CMD_ADDTAGPIG : m_lattice->addTaggedPairIG(); break;
    case CMD_GETNUMPARTICLES : getNumParticles();break;
    case CMD_ADDSIG : m_lattice->addSingleIG(); break;
    case CMD_ADDDAMP : m_lattice->addDamping(); break;
    case CMD_ADDTRIMESHIG : m_lattice->addTriMeshIG(); break;
    case CMD_ADDTRIMESH : m_lattice->addTriMesh(); break;
    case CMD_ADDMESH2D : m_lattice->addMesh2D(); break;
    case CMD_ADDMESH2DIG : m_lattice->addMesh2DIG(); break;
    case CMD_MAKELATTICE : makeLattice(); break;
    case CMD_INITLATTICE : initLattice(); break;
    case CMD_INITLATTICECIRC : initLatticeCirc(); break;
    case CMD_INITCOMPLEX : m_lattice->initComplex(); break;
    case CMD_ADDWALL : m_lattice->addWall(); break;
    case CMD_ADDEWALLIG : m_lattice->addElasticWIG(); break;
    case CMD_ADDBWALLIG : m_lattice->addBondedWIG(); break;
    case CMD_ADDVWALLIG : m_lattice->addViscWIG(); break;
    case CMD_ADDBBWALLIG : m_lattice->addDirBondedWIG(); break;
    case CMD_ADDTAGGEDEWALLIG : m_lattice->addTaggedElasticWIG(); break;
    case CMD_EXIG : m_lattice->setExIG(); break;
    case CMD_ADD_SPF : m_lattice->addScalarParticleField(); break;
    case CMD_ADD_VPF : m_lattice->addVectorParticleField(); break;
    case CMD_ADD_SIF : m_lattice->addScalarInteractionField(); break;
    case CMD_ADD_HIF : m_lattice->addScalarHistoryInteractionField(); break;
    case CMD_ADD_VIF : m_lattice->addVectorInteractionField(); break;
    case CMD_ADD_VTF : m_lattice->addVectorTriangleField(); break;
    case CMD_ADD_STF : m_lattice->addScalarTriangleField(); break;
    case CMD_ADD_VWF : m_lattice->addVectorWallField(); break;
    case CMD_SEND_FIELDS : m_lattice->sendFieldData(); break;
    case CMD_RECEIVEPARTICLES : m_lattice->receiveParticles();break;
    case CMD_RECEIVECONNECTIONS : m_lattice->receiveConnections();break;
    case CMD_ADDSPHEREBODY : m_lattice->addSphereBody();break;
    case CMD_ADDESPHEREBODYIG : m_lattice->addESphereBodyIG();break;
    case CMD_GETSPHEREBODYPOS : m_lattice->getSphereBodyPos();break;
    case CMD_GETSPHEREBODYFORCE : m_lattice->getSphereBodyForce();break;
    case CMD_SAVECHECKPOINT :
      {
        m_timersPtr->start("CheckPoint");
        m_pCheckPointer->saveRestartable();
        m_timersPtr->stop("CheckPoint");
        break;
      }
    case CMD_SAVESNAPSHOT :
      {
        m_timersPtr->start("SnapShot");
        m_pSnapShooter->saveDump();
        m_timersPtr->stop("SnapShot");
        break;
      }
    case CMD_SAVECHECKPOINTWTM :
      {
        m_timersPtr->start("CheckPointWTM");
        m_pCheckPointer->saveThroughMaster(m_tml_global_comm);
        m_timersPtr->stop("CheckPointWTM");
        break;
      }
    case CMD_LOADCHECKPOINT :
      {
        m_timersPtr->start("LoadCheckPoint");
        m_pCheckPointer->loadCheckPoint();
        m_timersPtr->stop("LoadCheckPoint");
        break;
      }
    case CMD_PTAG : m_lattice->tagParticleNearestTo();break;
    case CMD_FINDNEARESTPARTICLE: findParticleNearestToPoint();break;
    case CMD_GETPARTICLEPOSN: getParticlePosn();break;
    case CMD_GETWALLPOS: m_lattice->getWallPos();break;
    case CMD_GETWALLFORCE: m_lattice->getWallForce();break;
    case CMD_DO2DCALCULATIONS : do2dCalculations();break;
    case CMD_PERFORMTIMING : performTiming();break;
    case CMD_SAVETIMINGDATA : saveTimingData();break;
    case CMD_MOVENODE : m_lattice->moveSingleNode();break;
    case CMD_MOVETAGGEDNODES : m_lattice->moveTaggedNodes();break;
    case CMD_TRANSLATEMESHBY : translateMeshBy();break;
    case CMD_ROTATEMESHBY : rotateMeshBy();break;
    case CMD_ADDBONDEDTRIMESHIG : m_lattice->addBondedTriMeshIG(); break;
    case CMD_ADDBONDEDMESH2DIG : m_lattice->addBondedMesh2DIG() ; break;
    case CMD_REMOVEIG : m_lattice->removeIG() ; break;
    case CMD_GETMESHNODEREF : m_lattice->getMeshNodeRef() ; break;
    case CMD_GETMESHFACEREF : m_lattice->getMeshFaceRef() ; break;
    case CMD_GETMESH2DSTRESS : m_lattice->getMesh2DStress() ; break;
    case CMD_GETTRIMESHFORCE : m_lattice->getTriMeshForce() ; break;
    case CMD_GETIDPARTICLEDATA : getIdParticleData(); break;
    case CMD_IDPARTICLEMOVE : moveSingleParticle(); break;
    case CMD_SETTIMESTEPSIZE: setTimeStepSize(); break;
    case CMD_SETVERBOSITY : setVerbosity(); break;
	case CMD_INITCONSOLE : initializeConsole(); break;
	case CMD_SETCONSOLEFNAME : setConsoleFilename(); break;
	case CMD_SETCONSOLEBUFF : setConsoleBuffered(); break;
	case CMD_SETINTERACTIONPARAMS : m_lattice->setInteractionParameter(); break;
    case CMD_FINISH : break;
    default:
      {
        console.Error()
          << m_global_rank << " got unknown command : "
          << command << "\n";
        is_error=true;
      }
    }

    std::stringstream msg;
    msg << "cmd=" << command;
    m_tml_global_comm.barrier(msg.str().c_str());
    //m_tml_global_comm.barrier();
  } while ((command!=CMD_FINISH)&&(!is_error));
#ifdef __TIMING
  m_lattice->printTimes();
#endif //__TIMING
  m_timersPtr->clear();
  m_tml_global_comm.barrier("Final");
#ifdef HAVE_MPI_COMM_SPAWN
  MPI_Comm_disconnect(&m_global_comm);
#endif
}
