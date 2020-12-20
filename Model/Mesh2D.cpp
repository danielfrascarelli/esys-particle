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
#include "Foundation/console.h"
#include "Model/Mesh2D.h"
#include "Parallel/Mesh2DReader.h"

// --- STL includes ---
#include <set>

using std::set;
using namespace esys::lsm;
/*!
  constructor for empty 2D mesh
*/
Mesh2D::Mesh2D()
{}

/*!
  setup 2D mesh from node and edge data

  \param node_vec the node data
  \param edge_vec the edge data
*/ 
void Mesh2D::LoadMesh(const vector<MeshNodeData2D>& node_vec,const vector<MeshEdgeData2D>& edge_vec)
{
  // get temporary set of corner ids used in edges
  set<int> tmp_cid_set;
  for(vector<MeshEdgeData2D>::const_iterator iter=edge_vec.begin();
      iter!=edge_vec.end();
      iter++){
    int id0=iter->p1;
    int id1=iter->p2;
    tmp_cid_set.insert(id0);
    tmp_cid_set.insert(id1);
  }
  // setup node map
  for(vector<MeshNodeData2D>::const_iterator iter=node_vec.begin();
      iter!=node_vec.end();
      iter++){
    // check if corner is used in any edge
    set<int>::iterator sit=tmp_cid_set.find(iter->id);
    if(sit!=tmp_cid_set.end()){
      Corner2D Co=Corner2D(Vec3(iter->x,iter->y,0.0),iter->id);
      m_corners.push_back(Co);
      m_corner_by_id.insert(make_pair(iter->id,m_corners.size()-1));
    }
  }  
  console.Debug() << m_corner_by_id.size() << " elements in node map\n";
  // setup edges
  for(vector<MeshEdgeData2D>::const_iterator iter=edge_vec.begin();
      iter!=edge_vec.end();
      iter++){
    int id0=iter->p1;
    int id1=iter->p2;
    Vec3 p1=m_corners[m_corner_by_id[id0]].getPos();
    Vec3 p2=m_corners[m_corner_by_id[id1]].getPos();
    Edge2D E=Edge2D(id0,id1,p1,p2,iter->id,iter->tag);
    m_edges.push_back(E);
  }
  console.Debug() << m_edges.size() << " edges generated\n";
  // setup mapping from node id to edge
  for(vector<Edge2D>::iterator iter=m_edges.begin();
      iter!=m_edges.end();
      iter++){
    Edge2D* e_ptr=&(*iter);
    // get node ids
    int id0=(iter->getP0()).first;
    int id1=(iter->getP1()).first;
    // add pairs
    m_edge_by_node_id.insert(make_pair(id0,e_ptr));
    m_edge_by_node_id.insert(make_pair(id1,e_ptr));    
    // add edge pointer to corners
    m_corners[m_corner_by_id[id0]].addEdge(e_ptr);
    m_corners[m_corner_by_id[id1]].addEdge(e_ptr);
  }  
  // edge id to edge
  for(size_t i=0;i<m_edges.size();i++){
    Edge2D et=m_edges[i];
    //console.Debug() << et << "\n";
    m_edge_index_by_id[m_edges[i].getID()]=i;
  }

}

/*!
  move one node by a given amount

  \param id the id of the node
  \param d the displacement
*/
void Mesh2D::moveNode(int id,const Vec3& d)
{
  // move edges
  typedef multimap<int,Edge2D*>::iterator emmi;
  pair<emmi,emmi> efound=m_edge_by_node_id.equal_range(id);
  for(emmi iter=efound.first;iter!=efound.second;iter++){
    iter->second->moveNode(id,d);
  }
  // move nodes
  Corner2D *c=getCornerById(id);
  if(c!=NULL){
    c->move(d);
  }
}


void Mesh2D::translateBy(const Vec3 &translation)
{
  for (
    std::vector<Corner2D>::iterator it = m_corners.begin();
    it != m_corners.end();
    ++it
  )
  {
    it->move(translation);
  }
}


/*!
  zero all forces on the mesh. Currently forces are only on 
  edge
*/
void Mesh2D::zeroForces()
{
  for(vector<Edge2D>::iterator iter=m_edges.begin();
      iter!=m_edges.end();
      iter++){
    iter->zeroForce();
  }
}
  
/*!
  Get a pointer to a edge with a given ID. If the ID doesn't exist, return NULL

  \param id the id
*/
Edge2D* Mesh2D::getEdgeById(int id)
{
  Edge2D* res;
  map<int,int>::iterator it=m_edge_index_by_id.find(id);
  if(it!=m_edge_index_by_id.end()){
    res=&(m_edges[it->second]);
  } else {
    res=NULL;
  }

  return res;
}

/*!
  Get a pointer to a corner with a given ID. If the ID doesn't exist, return NULL

  \param id the id
*/
Corner2D* Mesh2D::getCornerById(int id)
{
  Corner2D* res;
  
  map<int,int>::iterator it=m_corner_by_id.find(id);
  if(it!=m_corner_by_id.end()){
    res=&(m_corners[it->second]);
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

  \warning doesn't deal with tags yet
*/
void Mesh2D::writeCheckPoint(ostream& ost,const string& delim) const
{
  int tag=0; 
  // build node map
  map<int,Vec3> nodemap;
  for(vector<Edge2D>::const_iterator iter=m_edges.begin();
      iter!=m_edges.end();
      iter++){
    nodemap.insert(iter->getP0());
    nodemap.insert(iter->getP1());    
  }
  // write nodes
  ost << "2D-Nodes " << nodemap.size() << delim;
  for(map<int,Vec3>::iterator iter=nodemap.begin();
      iter!=nodemap.end();
      iter++){
    ost << iter->first << " " << iter->first << " " << tag << " " 
	<< iter->second.X() << " " << iter->second.Y() << delim;
  }
  // write lines
  ost << "Line2 " << m_edges.size() << delim;
  int count=0;
  for(vector<Edge2D>::const_iterator iter=m_edges.begin();
      iter!=m_edges.end();
      iter++){
    ost << count << " " << tag << " ";
    ost << (iter->getP0()).first << " "; 
    ost << (iter->getP1()).first << delim;
    count++;
  }
}

/*!
  load checkpoint data from stream. Re-uses code from meshreader to read in 
  checkpointed meshes 

  \param ist the input stream

  \warning doesn't deal with tags yet
*/
void Mesh2D::loadCheckPoint(istream& ist)
{
  vector<MeshNodeData2D> node_vec;
  vector<MeshEdgeData2D> tri_vec;

  // --- Nodes ---
  Node2DReader nreader(ist);
  Node2DReader::Iterator &niter=nreader.getIterator();
  // read nodes into vector
  while(niter.hasNext()){
    node_vec.push_back(niter.next());
  }
 
  // --- Edges ---
  Edge2DReader treader(ist);
  Edge2DReader::Iterator &titer=treader.getIterator();
  // read edges into vector
  while(titer.hasNext()){
    tri_vec.push_back(titer.next());
  }
  LoadMesh(node_vec,tri_vec);
}
