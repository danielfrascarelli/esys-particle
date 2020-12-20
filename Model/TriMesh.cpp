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

#include "TriMesh.h"
#include "console.h"
#include "MeshReader.h"

using namespace esys::lsm;

/*!
  Constructor
*/
TriMesh::TriMesh()
{}

/*!
  setup the triangle mesh from node,edge and triangle data

  \param node_vec the node data
  \param tri_vec the triangle data
*/
void TriMesh::LoadMesh(const vector<MeshNodeData>& node_vec, const vector<MeshTriData>& tri_vec)
{

  console.XDebug() << " TriMesh::LoadMesh()\n";

  // setup node map
  for(vector<MeshNodeData>::const_iterator iter=node_vec.begin();
      iter!=node_vec.end();
      iter++){
    Corner Co=Corner(Vec3(iter->x,iter->y,iter->z),iter->id,iter->tag);
    m_corners.push_back(Co);
    m_corner_by_id.insert(make_pair(iter->id,m_corners.size()-1));
  }
  console.XDebug() << m_corner_by_id.size() << " elements in node map\n";

  // generate edges and connect them to the correct triangles
  map<pair<int,int>,vector<int> > edge_mmap;
  // generate triangles 
  for(vector<MeshTriData>::const_iterator iter=tri_vec.begin();
      iter!=tri_vec.end();
      iter++){
    int id0=iter->p1;
    int id1=iter->p2;
    int id2=iter->p3;
    Vec3 p1=m_corners[m_corner_by_id[id0]].getPos();
    Vec3 p2=m_corners[m_corner_by_id[id1]].getPos();
    Vec3 p3=m_corners[m_corner_by_id[id2]].getPos();
    Triangle Tri=Triangle(id0,id1,id2,p1,p2,p3,iter->id,iter->tag);
    m_triangles.push_back(Tri);
	
    // add edges to mmap
    int nt=m_triangles.size()-1;
    // 1-2 edge
    if(iter->p1<iter->p2){
      vector<int> &tv=(edge_mmap[make_pair(iter->p1,iter->p2)]);
      tv.push_back(nt);
    } else {
      vector<int> &tv=(edge_mmap[make_pair(iter->p2,iter->p1)]);
      tv.push_back(nt);
    }   
    // 2-3 edge
    if(iter->p2<iter->p3){
      vector<int> &tv=(edge_mmap[make_pair(iter->p2,iter->p3)]);
      tv.push_back(nt);
    } else {
      vector<int> &tv=(edge_mmap[make_pair(iter->p3,iter->p2)]);
      tv.push_back(nt);
    }  
    // 3-1 edge
    if(iter->p1<iter->p3){
      vector<int> &tv=(edge_mmap[make_pair(iter->p1,iter->p3)]);
      tv.push_back(nt);
    } else {
      vector<int> &tv=(edge_mmap[make_pair(iter->p3,iter->p1)]);
      tv.push_back(nt);
    }
  }
  
  // setup mapping from node id to triangle 
  for(vector<Triangle>::iterator iter=m_triangles.begin();
      iter!=m_triangles.end();
      iter++){
    Triangle* tr_ptr=&(*iter);
    // get node ids
    int id0=iter->getPid0();
    int id1=iter->getPid1();
    int id2=iter->getPid2();
    // add pairs
    m_triangle_by_node_id.insert(make_pair(id0,tr_ptr));
    m_triangle_by_node_id.insert(make_pair(id1,tr_ptr));
    m_triangle_by_node_id.insert(make_pair(id2,tr_ptr));
    // get indices of corners & add pointer to triangle
	m_corners[m_corner_by_id[id0]].addTriangle(tr_ptr);  
	m_corners[m_corner_by_id[id1]].addTriangle(tr_ptr);  
	m_corners[m_corner_by_id[id2]].addTriangle(tr_ptr);  
  }
  // and triangle id to triangle
  for(size_t i=0;i<m_triangles.size();i++){
    m_tri_index_by_id[m_triangles[i].getID()]=i;
  }
  console.XDebug() << m_triangles.size() << " triangles generated\n";

  // extract edges from multimap
  for(map<pair<int,int>,vector<int> >::iterator iter=edge_mmap.begin();
      iter!=edge_mmap.end();
      iter++){
    int id1=iter->first.first;
    int id2=iter->first.second;
    Vec3 p1=m_corners[m_corner_by_id[id1]].getPos();
    Vec3 p2=m_corners[m_corner_by_id[id2]].getPos();
    size_t ntri=(iter->second).size();
    switch(ntri){
    case 1: 
      {
	m_edges.push_back(Edge(id1,id2,p1,p2,&(m_triangles[iter->second[0]])));
      }; 
      break;
    case 2:
      {
	m_edges.push_back(Edge(id1,id2,p1,p2,&(m_triangles[iter->second[0]]),&(m_triangles[iter->second[1]])));
      }; 
      break;
    default : console.Error() << "Edge contained in " << ntri << " triangles ??"; break;
    }
  }
  // setup mapping from node id to edge
  for(vector<Edge>::iterator iter=m_edges.begin();
      iter!=m_edges.end();
      iter++){
    Edge* e_ptr=&(*iter);
    // get node ids
    pair<int,int> ids=iter->getIDs();
    // add pairs
    m_edge_by_node_id.insert(make_pair(ids.first,e_ptr));
    m_edge_by_node_id.insert(make_pair(ids.second,e_ptr));
  }

  console.XDebug() << "end TriMesh::LoadMesh()\n";
}

/*!
  Move a node in the mesh. If the node with the given Id isn't in the mesh,
  nothing happens

  \param id the id of the node
  \param d the displacement
*/
void TriMesh::moveNode(int id,const Vec3& d)
{
  // move triangles
  typedef multimap<int,Triangle*>::iterator tmmi;
  pair<tmmi,tmmi> tfound=m_triangle_by_node_id.equal_range(id);
  for(tmmi iter=tfound.first;iter!=tfound.second;iter++){
    iter->second->moveNode(id,d);
  }
  // move edges
  typedef multimap<int,Edge*>::iterator emmi;
  pair<emmi,emmi> efound=m_edge_by_node_id.equal_range(id);
  for(emmi iter=efound.first;iter!=efound.second;iter++){
    iter->second->moveNode(id,d);
  }
  // move corners
  m_corners[m_corner_by_id[id]].move(d);
}

void TriMesh::translateBy(const Vec3 &translation)
{
  // move triangles
  for(triangle_iterator iter=m_triangles.begin(); 
      iter!=m_triangles.end();
      iter++){
    iter->move(translation);
  }
  // move edges
  for(edge_iterator iter=m_edges.begin(); 
      iter!=m_edges.end();
      iter++){
    iter->move(translation);
  }
  // move corners
  for(corner_iterator it = m_corners.begin();
      it != m_corners.end();
      ++it){
    it->move(translation);
  }
}

/*!
    rotate the whol mesh

    \param origin a point on the rotation axis
    \param axis the orientation of the rotation axis
    \param angle the angle in radians
*/
void TriMesh::rotateBy(const Vec3& origin, const Vec3 &axis, double angle)
{
    // move triangles
    for(triangle_iterator iter=m_triangles.begin(); 
        iter!=m_triangles.end();
        iter++){
        iter->rotate(origin,axis,angle);
    }
    // move edges
    for(edge_iterator iter=m_edges.begin(); 
        iter!=m_edges.end();
        iter++){
        iter->rotate(origin,axis,angle);
    }
    // move corners
    for(corner_iterator iter = m_corners.begin();
        iter != m_corners.end();
        iter++){
        iter->rotate(origin,axis,angle);
    }
}

/*!
  check if any point in the mesh has moved by at least the given distance

  \param dist the distance 
*/
bool TriMesh::hasMovedBy(double dist)
{
  bool not_moved_by=true;

  std::vector<Corner>::iterator it = m_corners.begin();
  while ((it != m_corners.end()) && (not_moved_by)){
    // while (it != m_corners.end()){
    double dist_moved=it->getDistMoved();
    not_moved_by=(dist_moved<dist);
    it++;
  }
  
  return !not_moved_by;
}

/*!
  reset displacement since last neighbor search
*/
void TriMesh::resetCurrentDisplacement()
{
  for(std::vector<Corner>::iterator iter = m_corners.begin();
      iter!=m_corners.end();
      iter++){
    iter->resetOldPos();
  }
}



/*!
  Get a pointer to a triangle with a given ID. If the ID doesn't exist, return NULL

  \param id the id
*/
Triangle* TriMesh::getTriangleById(int id)
{
  Triangle* res;
  map<int,int>::iterator it=m_tri_index_by_id.find(id);
  if(it!=m_tri_index_by_id.end()){
    res=&(m_triangles[it->second]);
  } else {
    res=NULL;
  }

  return res;
}


/*!
  Write checkpoint data to stream. The mesh data is written in the original mesh file
  format -> can reuse meshreader to read in checkpointed meshes 

  \param ost the output stream
  \param delim the delimiter
*/
void TriMesh::writeCheckPoint(ostream& ost,const string& delim) const
{
  int tag=0; 
  // build node map
  map<int,Vec3> nodemap;
  for(vector<Triangle>::const_iterator iter=m_triangles.begin();
      iter!=m_triangles.end();
      iter++){
    nodemap.insert(iter->getP0());
    nodemap.insert(iter->getP1());
    nodemap.insert(iter->getP2());
  }
  // write nodes
  ost << "3D-Nodes " << nodemap.size() << delim;
  for(map<int,Vec3>::iterator iter=nodemap.begin();
      iter!=nodemap.end();
      iter++){
    ost << iter->first << " " << iter->first << " " << tag << " " << iter->second << delim;
  }
  ost << "Tri3 " << m_triangles.size() << delim;
  for(vector<Triangle>::const_iterator iter=m_triangles.begin();
      iter!=m_triangles.end();
      iter++){
    ost << iter->getID() << " " << iter->getTag() << " ";
    ost << (iter->getP0()).first << " "; 
    ost << (iter->getP1()).first << " "; 
    ost << (iter->getP2()).first << delim; 
  }
}

void TriMesh::loadCheckPoint(istream& ist)
{
  vector<MeshNodeData> node_vec;
  vector<MeshTriData> tri_vec;

  // --- Nodes ---
  NodeReader nreader(ist);
  NodeReader::Iterator &niter=nreader.getIterator();
  // read nodes into vector
  while(niter.hasNext()){
    node_vec.push_back(niter.next());
  }
 
  // --- Triangles ---
  TriReader treader(ist);
  TriReader::Iterator &titer=treader.getIterator();
  // read triangles into vector
  while(titer.hasNext()){
    tri_vec.push_back(titer.next());
  }
  LoadMesh(node_vec,tri_vec);
}

/*!
  zero all forces on the mesh. 
*/
void TriMesh::zeroForces()
{
  console.XDebug() << "TriMesh::zeroForces()\n";
  // triangles
  for(vector<Triangle>::iterator iter=m_triangles.begin();
      iter!=m_triangles.end();
      iter++){
    iter->zeroForce();
  }
  
  console.XDebug() << "end TriMesh::zeroForces()\n";
}

