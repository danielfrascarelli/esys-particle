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

// --- Project includes ---
#include "Parallel/GetRef_cmd.h"
#include "Parallel/GeometryReader.h"

// --- STL includes ---
#include <map>
#include <set>

using namespace esys::lsm;

/*!
  read model geometry from file

  \param fileName the name of the input file
*/
template <class TmplParticle>
void  CLatticeMaster::readGeometry(const std::string &fileName)
{
  console.Debug()
    << "begin CLatticeMaster::readGeometry: fileName="
    << fileName << "\n";

  // setup geometry reader and get geometry info from header file
  GeometryReader geoReader(fileName);
  GeometryInfo geoInfo = geoReader.getGeometryInfo();
  	
  // check if geometry info has been set previously
  if(m_bbx_has_been_set){ 
    // check if old & new geometry are compatible
    if(!m_geo_info.isCompatible(geoInfo)){
	std::cerr << "Geometry info read from file is incompatible with previously set geometry (bounding box, circular boundaries) - Model may not run properly!" << std::endl; 
    }
    // set only the geometry version, the other geometry info remains unchanged
    m_geo_info.setLsmGeoVersion(geoInfo.getLsmGeoVersion());
  } else {
    m_geo_info=geoInfo; // set full geometry meta-info

    // set model spatial domain -> sets geometry info & initializes neighbor table in workers
    if(m_geo_info.hasAnyPeriodicDimensions()){
	setSpatialDomain(m_geo_info.getBBoxCorners()[0],m_geo_info.getBBoxCorners()[1],m_geo_info.getPeriodicDimensions());
    } else {
	setSpatialDomain(m_geo_info.getBBoxCorners()[0],m_geo_info.getBBoxCorners()[1]);
    }
    m_bbx_has_been_set=true; // flag that spatial domain has been set
  }
  
  // read particles 
  GeometryReader::ParticleIterator &particleIt = geoReader.getParticleIterator();
  console.XDebug()<< "Number of particles: " << particleIt.getNumRemaining() << "\n";
  if (geoReader.getParticleType() != "Simple") {
    throw 
      std::runtime_error(
        (std::string("Unknown particle type ") + geoReader.getParticleType()).c_str()
      );
  }

  addParticles<GeometryReader::ParticleIterator,TmplParticle>(particleIt);

  // read connections
  GeometryReader::ConnectionIterator &connectionIt =
    geoReader.getConnectionIterator();
  console.XDebug()<< "Number of connections: " << connectionIt.getNumRemaining() << "\n";

  addConnections(connectionIt);
  console.Debug()<<"end CLatticeMaster::readGeometry\n";
}

template <typename TmplVisitor>
void CLatticeMaster::visitMeshFaceReferences(const string &meshName)
{
  throw std::runtime_error("CLatticeMaster::visitMeshFaceReferences: Not implemented.");
}

template <typename TmplVisitor>
void CLatticeMaster::visitMesh2dNodeReferences(const string &meshName, TmplVisitor &visitor)
{
  console.XDebug()<<"CLatticeMaster::visitMesh2dNodeReferences( " << meshName << ")\n";
  GetNodeRefCommand cmd(getGlobalRankAndComm(),meshName);
  // broadcast command
  cmd.broadcast();
  
  // receive data (multimap)
  std::multimap<int,int> ref_mmap;
  m_tml_global_comm.gather(ref_mmap);
  // collate into set
  std::set<int> ref_set; //== this is the set of node ids == 
  for(std::multimap<int,int>::iterator iter=ref_mmap.begin();
      iter!=ref_mmap.end();
      iter++){
    ref_set.insert(iter->second);
  }
  for (std::set<int>::const_iterator it = ref_set.begin(); it != ref_set.end(); it++)
  {
    visitor.visitNodeRef(*it);
  }
  console.XDebug()<<"end CLatticeMaster::visitMesh2dNodeReferences()\n";
}

template <typename TmplVisitor>
void CLatticeMaster::visitMesh2dEdgeStress(const string &meshName, TmplVisitor &visitor)
{
  console.XDebug()<<"CLatticeMaster::visitMesh2dEdgeStress( " << meshName << ")\n";
  std::multimap<int,pair<int,Vec3> > temp_mm;
  std::map<int,Vec3> data; //=== map of id, value ===

  BroadcastCommand cmd(getGlobalRankAndComm(), CMD_GETMESH2DSTRESS);
  cmd.append(meshName.c_str());
  cmd.broadcastCommand();
  cmd.broadcastBuffer();

  // get data from slaves
  m_tml_global_comm.gather(temp_mm);
  
  // add data together
  for(std::multimap<int,pair<int,Vec3> >::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    if(data.find((iter->second).first)==data.end()){ // id not in data -> insert
      data.insert(iter->second);
    } else { // id is in data -> add
      data[(iter->second).first]+=(iter->second).second;
    }
  }

  for (std::map<int,Vec3>::const_iterator it = data.begin(); it != data.end(); it++)
  {
    visitor.visitRefStressPair(it->first, it->second);
  }

  cmd.wait("visitMesh2dEdgeStress");
  console.XDebug()<<"end CLatticeMaster::visitMesh2dEdgeStress()\n";  
}

template <typename TmplVisitor>
void CLatticeMaster::visitTriMeshFaceForce(
  const string &meshName,
  TmplVisitor &visitor
)
{
  console.XDebug()<<"CLatticeMaster::visitTriMeshFaceForce( " << meshName << ")\n";
  std::multimap<int,pair<int,Vec3> > temp_mm;
  std::map<int,Vec3> data; //=== map of id, value ===

  BroadcastCommand cmd(getGlobalRankAndComm(), CMD_GETTRIMESHFORCE);
  cmd.append(meshName.c_str());
  cmd.broadcastCommand();
  cmd.broadcastBuffer();

  // get data from slaves
  m_tml_global_comm.gather(temp_mm);
  
  // add data together
  for(std::multimap<int,pair<int,Vec3> >::iterator iter=temp_mm.begin();
      iter!=temp_mm.end();
      iter++){
    if(data.find((iter->second).first)==data.end()){ // id not in data -> insert
      data.insert(iter->second);
    } else { // id is in data -> add
      data[(iter->second).first]+=(iter->second).second;
    }
  }

  for (std::map<int,Vec3>::const_iterator it = data.begin(); it != data.end(); it++)
  {
    visitor.visitRefForcePair(it->first, it->second);
  }

  cmd.wait("visitTriMeshFaceStress");
  console.XDebug()<<"end CLatticeMaster::visitTriMeshFaceForce()\n";  
}

template <typename TmplVisitor, typename TmplParticle>
void CLatticeMaster::visitParticlesOfType(
  const IdVector &particleIdVector,
  TmplVisitor &visitor
)
{
  console.Debug() << "CLatticeMaster::visitParticlesOfType: enter\n";
  typedef std::multimap<int,TmplParticle> ParticleMMap;
  ParticleMMap particleMMap;

  console.Debug()
    << "CLatticeMaster::visitParticlesOfType: broadcasting command\n";
  BroadcastCommand cmd(getGlobalRankAndComm(), CMD_GETIDPARTICLEDATA);
  cmd.broadcastCommand();

  console.Debug()
    << "CLatticeMaster::visitParticlesOfType: broadcasting particle id's\n";
  m_tml_global_comm.broadcast_cont(particleIdVector);

  console.Debug()
    << "CLatticeMaster::visitParticlesOfType:"
    << " gathering particle data from workers\n";
  m_tml_global_comm.gather_packed(particleMMap);
  console.Debug()
    << "CLatticeMaster::visitParticlesOfType:"
    << " gathered " << particleMMap.size() << " particles\n";

  console.Debug()
    << "CLatticeMaster::visitParticlesOfType:"
    << " visiting particle data\n";
  for(
    typename ParticleMMap::iterator iter=particleMMap.begin();
    iter != particleMMap.end();
    iter++
  )
  {
    iter->second.visit(visitor);
  }

  cmd.wait("visitParticles");
  console.Debug() << "CLatticeMaster::visitParticlesOfType: exit\n";
}

template <typename TmplVisitor>
void CLatticeMaster::visitParticles(
  const IdVector &particleIdVector,
  TmplVisitor &visitor
)
{
  console.Debug() << "CLatticeMaster::visitParticles: enter\n";

  if (m_particle_type == "Basic")
  {
    visitParticlesOfType<TmplVisitor,CParticle>(particleIdVector, visitor);
  }
  else if (m_particle_type == "Rot")
  {
    visitParticlesOfType<TmplVisitor,CRotParticle>(particleIdVector, visitor);
  }
  else if (m_particle_type == "RotVi")
  {
    visitParticlesOfType<TmplVisitor,CRotParticleVi>(particleIdVector, visitor);
  }
  else if (m_particle_type == "RotTherm")
  {
    visitParticlesOfType<TmplVisitor,CRotThermParticle>(particleIdVector, visitor);
  }
  else
  {
    throw
      std::runtime_error(
        std::string("Unknown particle type: ") + m_particle_type
      );
  }
  console.Debug() << "CLatticeMaster::visitParticles: exit\n";
}

template <class TmplIterator, class TmplParticle>
void CLatticeMaster::addParticles(TmplIterator &particleIt)
{
  CMPIBarrier barrier(m_global_comm);
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);

  // cast storage pointer to correct type 
  vector<TmplParticle> particleVector;
  console.XDebug()
    << "CLatticeMaster::addParticles:"
    << " Reserving vector memory for particles.\n";

  // Determine the directions for calculating the minimum and maximum extents of particles.
  // Indexes for valid directions x, y, z are 0, 1, 2.  If a dimension is irrelevant for the 
  // purpose of finding particle extents (because it is periodic or does not factor into 
  // a 2D calculation), its index is 3.  The indexes are sorted ascending to minimize the 
  // number of tests needed to process each particle in particlesMinMax(...).
  esys::lsm::IntVector periodicity(3);
  if (m_geo_info.hasAnyPeriodicDimensions())
    periodicity = m_geo_info.getPeriodicDimensions();
  for (int i=0; i<3; i++)
    m_particle_dimensions[i] = periodicity[i]==0 ? i : 3;
  if (m_geo_info.is2d())
    m_particle_dimensions[2]=3;
  int periodicTotal = periodicity[0] + periodicity[1] + periodicity[2];
  if (periodicTotal==1 || periodicTotal==2)
    std::sort(m_particle_dimensions.begin(), m_particle_dimensions.end());

  const int numBroadcastParticles = 50000;
  particleVector.reserve(numBroadcastParticles);
  console.XDebug()
    << "CLatticeMaster::addParticles:"
    << " Beginning add-particle loop..." << "\n";
  for (int i = 0; particleIt.hasNext(); i++) {
    const TmplParticle particle(particleIt.next());

    // Update the minimum and maximum extents of the particles as each is read in.
    particlesMinMax(particle);
    
    particleVector.push_back(particle);
    if (((i+1) % 5000) == 0) {
      console.XDebug()
        << "CLatticeMaster::addParticles:"
        << "Adding particle with id "
        << particleVector.rbegin()->getID() << "\n";
    }

    if (((i+1) % numBroadcastParticles) == 0)
    {
      console.XDebug() << "CLatticeMaster::addParticles:"
        << " Broadcasting receive cmd...." << "\n";
      cmd_buffer.broadcast(CMD_RECEIVEPARTICLES);
      console.XDebug()
        << "CLatticeMaster::addParticles: Broadcasting particles...." << "\n";
      m_tml_global_comm.broadcast_cont_packed(particleVector);
      barrier.wait("CLatticeMaster::addParticles: Post particle broadcast.");
      barrier.wait("CLatticeMaster::addParticles: Post Command.");
      particleVector.clear();
      particleVector.reserve(numBroadcastParticles);
    }
  }
  console.XDebug()
    << "CLatticeMaster::addParticles: Done add-particle loop..." << "\n";
  console.XDebug()
    << "CLatticeMaster::addParticles: Broadcasting final particle-receive cmd"
    << "\n";
  cmd_buffer.broadcast(CMD_RECEIVEPARTICLES);
  console.XDebug() << "CLatticeMaster::addParticles:"
    << " Broadcasting final set of particles...\n";
  m_tml_global_comm.broadcast_cont_packed(particleVector);
  barrier.wait("Post final particle broadcast.");
  barrier.wait("Post final particle-broadcast command.");
  // -- build ntable
  console.XDebug()
    << "CLatticeMaster::addParticles: "
    << "Building ntable (searchNeighbours)...\n";
  searchNeighbors(true);
  console.XDebug()
    << "CLatticeMaster::addParticles: exit\n";

}

template <class TmplIterator>
void CLatticeMaster::addConnections(TmplIterator &connectionIt)
{
  console.XDebug()
    << "CLatticeMaster::addConnections: enter\n";
  
  CMPIBarrier barrier(m_global_comm);
  CMPILCmdBuffer cmd_buffer(m_global_comm,m_global_rank);
  
  const int numBroadcastConnections = 100000;
  vector<int> connectionBuffer;
  connectionBuffer.reserve(numBroadcastConnections);
  int i = 0;
  for (; connectionIt.hasNext(); i++)
  {
    typename TmplIterator::value_type data = connectionIt.next();
    connectionBuffer.push_back(data.getTag());
    connectionBuffer.push_back(data.getP1Id());
    connectionBuffer.push_back(data.getP2Id());
    if ((i+1) % 50000 == 0)
    {
      console.XDebug() << "Adding connection number " << i << "\n";
    }
    if ((i+1) % numBroadcastConnections == 0)
    {
      console.XDebug() << "CLatticeMaster::addConnections:"
        << " Broadcasting receive cmd...." << "\n";
      cmd_buffer.broadcast(CMD_RECEIVECONNECTIONS);
      console.XDebug()
        << "CLatticeMaster::addConnections: Broadcasting connections...." << "\n";
      m_tml_global_comm.broadcast_cont_packed(connectionBuffer);
      barrier.wait("CLatticeMaster::addConnections: Post connection broadcast.");
      barrier.wait("CLatticeMaster::addConnections: Post Command.");
      connectionBuffer.clear();
      connectionBuffer.reserve(numBroadcastConnections);
    }
    //m_temp_conn[data.getTag()].push_back(data.getP1Id());
    //m_temp_conn[data.getTag()].push_back(data.getP2Id());
  }
  console.XDebug()
    << "CLatticeMaster::addConnections: Done add-connection loop..." << "\n";
  console.XDebug()
    << "CLatticeMaster::addConnections: Broadcasting final connection-receive cmd"
    << "\n";
  cmd_buffer.broadcast(CMD_RECEIVECONNECTIONS);
  console.XDebug() << "CLatticeMaster::addConnections:"
    << " Broadcasting final set of connections...\n";
  m_tml_global_comm.broadcast_cont_packed(connectionBuffer);
  barrier.wait("Post final connection broadcast.");
  barrier.wait("Post final connection-broadcast command.");

  console.XDebug()<< "Added " << i << " connections..." << "\n";
  // -- build ntable
  searchNeighbors(true);
  console.XDebug()
    << "CLatticeMaster::addConnections: exit\n";
}

template<typename TmplParticle>
void CLatticeMaster::particlesMinMax(const TmplParticle &particle)
{
  Vec3 position = particle.getPos();
  double radius = particle.getRad();

  for (int j=0, k=m_particle_dimensions[j]; j<3 && k<3; k=m_particle_dimensions[++j]) {
    if (m_init_min_pt[k]!=m_init_min_pt[k] || position[k]-radius<m_init_min_pt[k])
      m_init_min_pt[k] = position[k]-radius;
    if (m_init_max_pt[k]!=m_init_max_pt[k] || position[k]+radius>m_init_max_pt[k])
      m_init_max_pt[k] = position[k]+radius;
  }
}
