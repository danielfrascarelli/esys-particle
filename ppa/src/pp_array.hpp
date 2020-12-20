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

///--- IO includes ---
#include "Foundation/console.h"
#include <iostream>
#include <sstream>
#include <stdexcept>

#ifdef __INTEL_COMPILER
#include <mathimf.h>
#else
#include <math.h>
#endif
using std::cout;
using std::cerr;
using std::endl;
using std::flush;

//--- TML includes ---
#include "ntable/src/nt_slab.h"

//--- STL includes ---
#include <utility>

using std::pair;
using std::make_pair;



template<typename T>
const int ParallelParticleArray<T>::m_exchg_tag=42;

/*!
  Construct a parallel particle array from a given communicator and geometry,
  i.e. minimum and maximum corners and search range. The process topology is generated
  from the communicator (via MPI_Dims_create). All boundaries are assumed to be open
  (i.e. not circular)

  \param comm the communicator
  \param dims the dimensions the process space.
              {dims[0]=0; dims[1]=0; dims[2]=0;}--is 3D allocation of processors.
	      {dims[0]=0; dims[1]=0; dims[2]=1;}--is 2D allocation of processors.
	      {dims[0]=0; dims[1]=1; dims[2]=1;}--is 1D allocation of processors.
  \param min the (global) minimum corner of the model space
  \param max the (global) maximum corner of the model space
  \param range the search range
  \param alpha the pair search cutoff

  \todo check for (dims.size() > 3)
*/
template<typename T>
ParallelParticleArray<T>::ParallelParticleArray(TML_Comm *comm, const std::vector<unsigned int> &dims,const Vec3& min, const Vec3& max ,double range, double alpha):AParallelParticleArray(comm, dims), m_nt(NULL)
{
  //- set x-edges to non-circular
  m_circ_edge_x_down=false;
  m_circ_edge_x_up=false;
  m_fluidinitiated=false; //fluid contents
  m_firststep=true; //fluid contents
  //- get own process coords
  vector<int> pcoords=m_comm.get_coords();
  //- initialize ntable
  // get global number of cells
  int nx_global,ny_global,nz_global;
  nx_global=lrint((max[0]-min[0])/range);
  ny_global=lrint((max[1]-min[1])/range);
  nz_global=lrint((max[2]-min[2])/range);
  if ((fabs((double(nx_global)-(max[0]-min[0])/range)) > 1e-6) ||
      (fabs((double(ny_global)-(max[1]-min[1])/range)) > 1e-6) ||
      (fabs((double(nz_global)-(max[2]-min[2])/range)) > 1e-6)){
    m_nt=NULL;
    std::stringstream msg;
    msg << "size doesn't fit range" << endl; // replace by throw
    msg << "diff x : " << double(nx_global)-(max[0]-min[0])/range << endl;
    msg << "diff y : " << double(ny_global)-(max[1]-min[1])/range << endl;
    msg << "diff z : " << double(nz_global)-(max[2]-min[2])/range << endl;
    console.Error() << msg.str() << "\n";
    throw std::runtime_error(msg.str());
  } else {
    // calc local min. cell, considering overlap
    int nx_min=((nx_global*pcoords[0])/m_comm.get_dim(0))-1;
    int ny_min=((ny_global*pcoords[1])/m_comm.get_dim(1))-1;
    int nz_min=((nz_global*pcoords[2])/m_comm.get_dim(2))-1;
    // calc local number of cells, considering overlap
    int nx=(((nx_global*(pcoords[0]+1))/m_comm.get_dim(0))-nx_min)+1;
    int ny=(((ny_global*(pcoords[1]+1))/m_comm.get_dim(1))-ny_min)+1;
    int nz=(((nz_global*(pcoords[2]+1))/m_comm.get_dim(2))-nz_min)+1;

    // init
    m_nt=new NeighborTable<T>(nx,ny,nz,range,alpha,min,max,nx_min,ny_min,nz_min);
  }
  // init time stamp
  m_timestamp=0;
}

/*!
  Construct a parallel particle array from a given communicator and geometr,
  i.e. minimum and maximum corners and search range. The process topology is generated
  from the communicator (via MPI_Dims_create). The boundary conditions i.e. circular or open
  are given as parameter.

  \param comm the communicator
  \param dims the dimensions the process space.
              {dims[0]=0; dims[1]=0; dims[2]=0;}--is 3D allocation of processors.
	      {dims[0]=0; dims[1]=0; dims[2]=1;}--is 2D allocation of processors.
	      {dims[0]=0; dims[1]=1; dims[2]=1;}--is 1D allocation of processors.
  \param circ circular/open boundary conditions
              {circ[0] : x-direction, circ[1] : y-direction and circ[2] : z-direction,
	      {true : circular, false : open}
  \param min the (global) minimum corner of the model space
  \param max the (global) maximum corner of the model space
  \param range the search range
  \param alpha the pair search cutoff

  \todo check for (dims.size() > 3)
*/
template<typename T>
ParallelParticleArray<T>::ParallelParticleArray(TML_Comm *comm, const vector<unsigned int> &dims, const vector<bool> &circ,const Vec3 &min,const Vec3 &max, double range, double alpha):AParallelParticleArray(comm,dims,circ), m_nt(NULL)
{
  m_fluidinitiated=false; //fluid contents
  m_firststep=true; //fluid contents
  //- get own process coords
  vector<int> pcoords=m_comm.get_coords();

  //- initialize ntable
  // get global number of cells
  int nx_global,ny_global,nz_global;
  nx_global=lrint((max[0]-min[0])/range);
  ny_global=lrint((max[1]-min[1])/range);
  nz_global=lrint((max[2]-min[2])/range);
  if ((fabs((double(nx_global)-(max[0]-min[0])/range)) > 1e-6) ||
      (fabs((double(ny_global)-(max[1]-min[1])/range)) > 1e-6) ||
      (fabs((double(nz_global)-(max[2]-min[2])/range)) > 1e-6)){
    m_nt=NULL;
    std::stringstream msg;
    msg << "size doesn't fit range" << endl; // replace by throw
    msg << "diff x : " << double(nx_global)-(max[0]-min[0])/range << endl;
    msg << "diff y : " << double(ny_global)-(max[1]-min[1])/range << endl;
    msg << "diff z : " << double(nz_global)-(max[2]-min[2])/range << endl;
    throw std::runtime_error(msg.str());
  } else {
    // calc local min. cell, considering overlap
    int nx_min=((nx_global*pcoords[0])/m_comm.get_dim(0))-1;
    int ny_min=((ny_global*pcoords[1])/m_comm.get_dim(1))-1;
    int nz_min=((nz_global*pcoords[2])/m_comm.get_dim(2))-1;
    // calc local number of cells, considering overlap
    int nx=(((nx_global*(pcoords[0]+1))/m_comm.get_dim(0))-nx_min)+1;
    int ny=(((ny_global*(pcoords[1]+1))/m_comm.get_dim(1))-ny_min)+1;
    int nz=(((nz_global*(pcoords[2]+1))/m_comm.get_dim(2))-nz_min)+1;
    // init local neighbortable
    m_nt=new NeighborTable<T>(nx,ny,nz,range,alpha,min,max,nx_min,ny_min,nz_min);
    // setup circular shift values
    m_xshift=(max-min).X();
    m_yshift=(max-min).Y();
    m_zshift=(max-min).Z();
    // setup circular edges
    m_circ_edge_x_down=(circ[0] && (pcoords[0]==0)) ? true : false ;
    m_circ_edge_x_up=(circ[0] && (pcoords[0]==m_comm.get_dim(0)-1)) ? true : false ;
  }
  // init time stamp
  m_timestamp=0;
}

/*!
  destructor
*/
template<typename T>
ParallelParticleArray<T>::~ParallelParticleArray()
{
  if(m_nt!=NULL) delete m_nt;
}


/**** fluid contents: begin ****/
template<typename T>
void ParallelParticleArray<T>::setFluid(double cellside,double Bw,double Bp,double Mu,double alpha,double flowrate,double pressure,Vec3 inflow,Vec3 outflow)
{
  Vec3 min=m_nt->base_point();
  Vec3 max=m_nt->max_point();
  double dim=m_nt->dim();
  int dimension;
  vector<int> pcoords=m_comm.get_coords();

  double xsize=max.X()-min.X();
  int nt_x=lrint(xsize/dim);
  dimension=m_comm.get_dim(0);
  int nx=lrint(nt_x/dimension*dim/cellside)+2;

  double ysize=max.Y()-min.Y();
  int nt_y=lrint(ysize/dim);
  dimension=m_comm.get_dim(1);
  int ny=lrint(nt_y/dimension*dim/cellside)+2;

  double zsize=max.Z()-min.Z();
  int nt_z=lrint(zsize/dim);
  dimension=m_comm.get_dim(2);
  int nz=lrint(nt_z/dimension*dim/cellside)+2;

  int nx_min=pcoords[0]*(nx-2)-1;
  int ny_min=pcoords[1]*(ny-2)-1;
  int nz_min=pcoords[2]*(nz-2)-1;

  m_nt->initiateFluid(Bw,Bp,Mu,alpha,flowrate,pressure,inflow,outflow,nx,ny,nz,nx_min,ny_min,nz_min);
  m_fluidinitiated=true;
}


template<typename T>
void ParallelParticleArray<T>::setFluidVec3(Vec3 cellside,double Bw,double Bp,double Mu,double alpha,double flowrate,double pressure,Vec3 inflow,Vec3 outflow)
{
  Vec3 min=m_nt->base_point();
  Vec3 max=m_nt->max_point();
  double dim=m_nt->dim();
  int dimension;
  vector<int> pcoords=m_comm.get_coords();

  double xsize=max.X()-min.X();
  int nt_x=lrint(xsize/dim);
  dimension=m_comm.get_dim(0);
  int nx=lrint(nt_x/dimension*dim/cellside[0])+2;

  double ysize=max.Y()-min.Y();
  int nt_y=lrint(ysize/dim);
  dimension=m_comm.get_dim(1);
  int ny=lrint(nt_y/dimension*dim/cellside[1])+2;

  double zsize=max.Z()-min.Z();
  int nt_z=lrint(zsize/dim);
  dimension=m_comm.get_dim(2);
  int nz=lrint(nt_z/dimension*dim/cellside[2])+2;

  int nx_min=pcoords[0]*(nx-2)-1;
  int ny_min=pcoords[1]*(ny-2)-1;
  int nz_min=pcoords[2]*(nz-2)-1;

  m_nt->initiateFluid(Bw,Bp,Mu,alpha,flowrate,pressure,inflow,outflow,nx,ny,nz,nx_min,ny_min,nz_min);
  m_fluidinitiated=true;
}


template<typename T>
void ParallelParticleArray<T>::updateFluid()
{
  if(m_firststep){
    m_nt->calPorositySphere0();
    m_firststep=false;
  }
  m_nt->calPorositySphere();
  m_nt->updateFluidcells();
}


template<typename T>
void ParallelParticleArray<T>::setPressure(vector<pair<Vec3,double> > pressure)
{
  m_nt->setPressure(pressure);
}


template<typename T>
void ParallelParticleArray<T>::exchangeCells(double timestep,int nt)
{
  // x-dir (there is at least one dimension)
  if(m_comm.get_dim(0)>1){
    // upstream
    exchange_boundary(m_nt->get_yz_cellslab(m_nt->xcell()-2),0,1);
    // downstream
    exchange_boundary(m_nt->get_yz_cellslab(1),0,-1);
  }

  // y-dir
  if(m_comm.get_dim(1)>1){
    // upstream
    exchange_boundary(m_nt->get_xz_cellslab(m_nt->ycell()-2),1,1);
    // downstream
    exchange_boundary(m_nt->get_xz_cellslab(1),1,-1);
  }

  // z-dir
  if(m_comm.get_dim(2)>1){
    // upstream
    exchange_boundary(m_nt->get_xy_cellslab(m_nt->zcell()-2),2,1);
    // downstream
    exchange_boundary(m_nt->get_xy_cellslab(1),2,-1);
  }

  m_nt->calCoeffi(timestep,nt);
  //m_nt->calVelocity();
}


template<typename T>
void ParallelParticleArray<T>::exchange_boundary(vector<CFluidCell> send_data, int dir, int dist)
{
  vector<CFluidCell> recv_data;
  // exchange
  m_comm.shift_cont_packed(send_data,recv_data,dir,dist,m_exchg_tag);
  // apply received data
  if(dir==0){
    if(dist==1){m_nt->set_yz_cellslab(0,recv_data);}
    else if(dist==-1){m_nt->set_yz_cellslab(m_nt->xcell()-1,recv_data);}
  }
  else if(dir==1){
    if(dist==1){m_nt->set_xz_cellslab(0,recv_data);}
    else if(dist==-1){m_nt->set_xz_cellslab(m_nt->ycell()-1,recv_data);}
  }
  else if(dir==2){
    if(dist==1){m_nt->set_xy_cellslab(0,recv_data);}
    else if(dist==-1){m_nt->set_xy_cellslab(m_nt->zcell()-1,recv_data);}
  }
}


/*!
  Get a value for each fluid cell using a fluid cell member function and return a vector
  of values, with global position of cell.

  \param rdf the fluid cell member function
*/
template<typename T>
template<typename P>
vector<pair<Vec3,P> > ParallelParticleArray<T>::forAllInnerCellsGet(P (CFluidCell::*rdf)() const)
{
  vector<pair<Vec3,P> > res;
  res=m_nt->forAllInnerCellsGet(rdf);
  return res;
}


/*!
  Get a value for each fluid cell using a fluid cell member function and return a vector
  of values, with global index of cell.

  \param rdf the fluid cell member function
*/
template<typename T>
template<typename P>
vector<pair<Vec3,P> > ParallelParticleArray<T>::forAllInnerCellsGetIndexed(P (CFluidCell::*rdf)() const)
{
  vector<pair<Vec3,P> > res;
  res=m_nt->forAllInnerCellsGetIndexed(rdf);
  return res;
}


template<typename T>
template<typename P>
vector<P> ParallelParticleArray<T>::forAllInnerCellsGetSum(P (CFluidCell::*rdf)() const)
{
  vector<P> res;
  res=m_nt->forAllInnerCellsGetSum(rdf);
  return res;
}

/**** fluid contents: end ****/


/*!
  insert a single particle into the storage

  \param p the particle
*/
template<typename T>
void ParallelParticleArray<T>::insert(const T& p)
{
  m_nt->insert(p);
}

/*!
  insert a STL vector of particles into the storage

  \param vp the vector of particles
*/
template<typename T>
void ParallelParticleArray<T>::insert(const vector<T>& vp)
{
  for(typename vector<T>::const_iterator iter=vp.begin();
      iter!=vp.end();
      iter++){
    m_nt->insert(*iter);
  }
}

/*!
  check if a position is in the inner part

  \param pos the position
*/
template<typename T>
bool ParallelParticleArray<T>::isInInner(const Vec3& pos)
{
  return m_nt->isInInner(pos);
}

/*!
  Get the pointer to a particle with a given id. Return NULL if
  there is no particle with this index.

  \param id the particle id.
*/
template<typename T>
T* ParallelParticleArray<T>::getParticlePtrByIndex(int id)
{
  return m_nt->ptr_by_id(id);
}

/*!
  Get the pointer to a particle closest to a given position. Return NULL if
  the position is outside the area.

  \param pos the position.
*/
template<typename T>
T* ParallelParticleArray<T>::getParticlePtrByPosition(const Vec3& pos)
{
  return m_nt->getNearestPtr(pos);
}
/*!
  Rebuild the neighbor table, i.e. relocate particles to the
  appropriate gridpoints and exchange boundary particles with
  neighboring nodes. No (geometric) resizing done.
*/
template<typename T>
void ParallelParticleArray<T>::rebuild()
{
  // cout << "PPA at node " << m_comm.rank() << "rebuilding " << endl << flush;
  //- get own process coords (for debug output)
  vector<int> pcoords=m_comm.get_coords();

  // rebuild locally
  // cout << "node " << m_comm.rank() << " pre-build " << *m_nt << endl;
  m_nt->build();
  // cout << "node " << m_comm.rank() << " pre-exchange " << *m_nt << endl;
  //- exchange boundary particles
  NTSlab<T> up_slab;
  NTSlab<T> down_slab;
  vector<T> recv_data;
  bool circ_edge=false;

  // x-dir (there is at least one dimension)
  if(m_comm.get_dim(0)>1){
    // clean out boundary slabs
    NTSlab<T> up_boundary_slab=m_nt->yz_slab(m_nt->xsize()-1);
    up_boundary_slab.erase(up_boundary_slab.begin(),up_boundary_slab.end());
    NTSlab<T> down_boundary_slab=m_nt->yz_slab(0);
    down_boundary_slab.erase(down_boundary_slab.begin(),down_boundary_slab.end());

    // synchronize
    m_comm.barrier();

    // define tranfer slabs
    up_slab=m_nt->yz_slab(m_nt->xsize()-2);
    down_slab=m_nt->yz_slab(1);

    // upstream
    if(m_circ_edge_x_up){ // circ. bdry
      // cout << "circular shift, node " << m_comm.rank() << " x-dim, up slab size " << up_slab.size() << endl;
      // copy particles from transfer slab into buffer
      vector<T> buffer;
      for(typename NTSlab<T>::iterator iter=up_slab.begin();
	  iter!=up_slab.end();
	  iter++){
	buffer.push_back(*iter);
      }
      // shift particles in buffer by circular shift value
      for(typename vector<T>::iterator iter=buffer.begin();
	  iter!=buffer.end();
	  iter++){
	iter->setCircular(Vec3(-1.0*m_xshift,0.0,0.0));
      }
      m_comm.shift_cont_packed(buffer,*m_nt,0,1,m_exchg_tag);
    } else {
      m_comm.shift_cont_packed(up_slab,*m_nt,0,1,m_exchg_tag);
    }
    // downstream
    if(m_circ_edge_x_down){// circ. bdry
      // cout << "circular shift , node " << m_comm.rank() << " x-dim, down slab size " << down_slab.size() << endl;
      // copy particles from transfer slab into buffer
      vector<T> buffer;
      for(typename NTSlab<T>::iterator iter=down_slab.begin();
	  iter!=down_slab.end();
	  iter++){
	buffer.push_back(*iter);
      }
      // shift particles in buffer by circular shift value
      for(typename vector<T>::iterator iter=buffer.begin();
	  iter!=buffer.end();
	  iter++){
	iter->setCircular(Vec3(m_xshift,0.0,0.0));
      }
      m_comm.shift_cont_packed(buffer,*m_nt,0,-1,m_exchg_tag);
    } else {
      m_comm.shift_cont_packed(down_slab,*m_nt,0,-1,m_exchg_tag);
    }
  }
  // y-dir
  if(m_comm.get_dim(1)>1){
    // clean out boundary slabs
    NTSlab<T> up_boundary_slab=m_nt->xz_slab(m_nt->ysize()-1);
    up_boundary_slab.erase(up_boundary_slab.begin(),up_boundary_slab.end());
    NTSlab<T> down_boundary_slab=m_nt->xz_slab(0);
    down_boundary_slab.erase(down_boundary_slab.begin(),down_boundary_slab.end());

    // synchronize
    m_comm.barrier();

    // define tranfer slabs
    up_slab=m_nt->xz_slab(m_nt->ysize()-2);
    down_slab=m_nt->xz_slab(1);
    if(!circ_edge){ // normal boundary
      // upstream
      m_comm.shift_cont_packed(up_slab,*m_nt,1,1,m_exchg_tag);
      // downstream
      m_comm.shift_cont_packed(down_slab,*m_nt,1,-1,m_exchg_tag);
    } else { // circular edge -> buffer & shift received particles
         cout << "circular y shift not implemented" << endl;
    }
  }
  // z-dir
  if(m_comm.get_dim(2)>1){
    // cout << "zdim= " << m_comm.get_dim(2) << " shifting" << endl;
    // clean out boundary slabs
    NTSlab<T> up_boundary_slab=m_nt->xy_slab(m_nt->zsize()-1);
    up_boundary_slab.erase(up_boundary_slab.begin(),up_boundary_slab.end());
    NTSlab<T> down_boundary_slab=m_nt->xy_slab(0);
    down_boundary_slab.erase(down_boundary_slab.begin(),down_boundary_slab.end());

    // define tranfer slabs
    up_slab=m_nt->xy_slab(m_nt->zsize()-2);
    down_slab=m_nt->xy_slab(1);
    if(!circ_edge){ // normal boundary
      // upstream
      m_comm.shift_cont_packed(up_slab,*m_nt,2,1,m_exchg_tag);
      // downstream
      m_comm.shift_cont_packed(down_slab,*m_nt,2,-1,m_exchg_tag);
    } else { // circular edge -> buffer & shift received particles
         cout << "circular x shift not implemented" << endl;
    }
  }
  // update timestamp
  m_timestamp++;
  //cout<<"  ppa rebuilt: time_stamp="<<m_timestamp<<endl;
  // debug
 //  cout << "node " << m_comm.rank() << "post-exchange " << *m_nt << flush;
}

/*!
  For all particles shared with neighboring nodes, exchange some
  value accessible by read and write functions.

  \param rdf the particle member function to read the value
  \param wrtf the particle member function to write the value
*/
template<typename T>
template<typename P>
void ParallelParticleArray<T>::exchange(P (T::*rdf)(),void (T::*wrtf)(const P&))
{
  // x-dir (there is at least one dimension)
  if(m_comm.get_dim(0)>1){
    // cout << "xdim= " << m_comm.get_dim(0) << " exchanging" << endl;
    // do transfer
    exchange_single(rdf,wrtf,m_nt->yz_slab(m_nt->xsize()-2),m_nt->yz_slab(0),0,1);
    // downstream
    exchange_single(rdf,wrtf,m_nt->yz_slab(1),m_nt->yz_slab(m_nt->xsize()-1),0,-1);
  }
  // y-dir
  if(m_comm.get_dim(1)>1){
    // cout << "ydim= " << m_comm.get_dim(1) << " shifting" << endl;
    // upstream
    exchange_single(rdf,wrtf,m_nt->xz_slab(m_nt->ysize()-2),m_nt->xz_slab(0),1,1);
    // downstream
    exchange_single(rdf,wrtf,m_nt->xz_slab(1),m_nt->xz_slab(m_nt->ysize()-1),1,-1);
  }
  // z-dir
  if(m_comm.get_dim(2)>1){
    // cout << "zdim= " << m_comm.get_dim(2) << " shifting" << endl;
    // upstream
    exchange_single(rdf,wrtf,m_nt->xy_slab(m_nt->zsize()-2),m_nt->xy_slab(0),2,1);
    // downstream
    exchange_single(rdf,wrtf,m_nt->xy_slab(1),m_nt->xy_slab(m_nt->zsize()-1),2,-1);
  }
}

/*!
  Helper function which does the actual shifting of values for exchange

  \param rdf  the particle member function to read the value
  \param wrtf the particle member function to write the value
  \param send_slab the NTSlab in which the data to be sent are located
  \param recv_slab the NTSlab in which the received data will be applied
  \param dir the direction of the transfer (x,y,z)
  \param dist the shift distance, i.e. up or down (1,-1)
*/
template<typename T>
template<typename P>
void ParallelParticleArray<T>::exchange_single(P (T::*rdf)(),void (T::*wrtf)(const P&),
					       NTSlab<T> send_slab,NTSlab<T> recv_slab,
					       int dir,int dist)
{
  vector<P> send_data;
  vector<P> recv_data;

  // get data
  for(typename NTSlab<T>::iterator iter=send_slab.begin();
      iter!=send_slab.end();
      iter++){
    send_data.push_back(((*iter).*rdf)());
//     cout << "put exchange values from particle " << iter->getID() << " into buffer" << endl;
  }
  // exchange
  m_comm.shift_cont_packed(send_data,recv_data,dir,dist,m_exchg_tag);

  // apply received data
  // check if sizes are correct
  if(recv_data.size()==recv_slab.size()){
    int count=0;
    for(typename NTSlab<T>::iterator iter=recv_slab.begin();
	iter!=recv_slab.end();
	iter++){
      ((*iter).*wrtf)(recv_data[count]);
//       cout << "wrote exchange value to particle " << iter->getID() << endl;
      count++;
    }
  }else{
    cerr << "rank = " << m_comm.rank() << ":ParallelParticleArray<T>::exchange_single ERROR: size mismatch in received data, recv_data.size()!=recv_slab.size() - (" << recv_data.size() << "!=" << recv_slab.size() << "). dir = " << dir << ", dist = " << dist << endl;
  }
}

/*!
  Call a member function taking no argument for one particle. Do nothing if the
  particle with the id is not in the ntable.

  \param id the id of the particle
  \param mf the member function

  \warning current implementation is O(n)
*/
template<typename T>
void ParallelParticleArray<T>::forParticle(int id,void (T::*mf)())
{
  T* pp=m_nt->ptr_by_id(id);
  if(pp!=NULL){
    (pp->*mf)();
  }
}

/*!
  Call a member function taking one argument for one particle. Do nothing if the
  particle with the id is not in the ntable.

  \param id the id of the particle
  \param mf the member function
  \param arg the argument to the function call

  \warning current implementation is O(n)
*/
template<typename T>
template<typename P>
void ParallelParticleArray<T>::forParticle(int id,void (T::*mf)(P),const P& arg)
{
  T* pp=m_nt->ptr_by_id(id);
  if(pp!=NULL){
    (pp->*mf)(arg);
  }
}

/*!
  Call a member function taking no argument for all particles with a given tag.

  \param tag the tag
  \param mf the member function
*/
template<typename T>
void ParallelParticleArray<T>::forParticleTag(int tag,void (T::*mf)())
{
  for(typename NeighborTable<T>::iterator iter=m_nt->begin();
      iter!=m_nt->end();
      iter++){
    if(iter->getTag()==tag){
      ((*iter).*mf)();
    }
  }
}

/*!
  Call a member function taking one argument for all particleswith a given tag.

  \param tag the tag
  \param mf the member function
  \param arg the argument to the function call
*/
template<typename T>
template<typename P>
void ParallelParticleArray<T>::forParticleTag(int tag,void (T::*mf)(P),const P& arg)
{
  for(typename NeighborTable<T>::iterator iter=m_nt->begin();
      iter!=m_nt->end();
      iter++){
    if(iter->getTag()==tag){
      ((*iter).*mf)(arg);
    }
  }
}

/*!
  Call a member function taking no argument for all particles with a given tag and mask.
  The functions is called if the masked bits in the particle tag and the given tag are identical,
  i.e. if ptag & mask == tag & mask

  \param tag the tag
  \param mask the mask
  \param mf the member function
*/
template<typename T>
void ParallelParticleArray<T>::forParticleTagMask(int tag,int mask, void (T::*mf)())
{
  for(typename NeighborTable<T>::iterator iter=m_nt->begin();
      iter!=m_nt->end();
      iter++){
    if((iter->getTag() & mask) == (tag & mask)){
      ((*iter).*mf)();
    }
  }
}

/*!
  Call a member function taking one argument for all particles with a given tag and mask.
  The functions is called if the masked bits in the particle tag and the given tag are identical,
  i.e. if ptag & mask == tag & mask

  \param tag the tag
  \param mask the mask
  \param mf the member function
  \param arg the argument to the function call
*/
template<typename T>
template<typename P>
void ParallelParticleArray<T>::forParticleTagMask(int tag,int mask,void (T::*mf)(P),const P& arg)
{
  for(typename NeighborTable<T>::iterator iter=m_nt->begin();
      iter!=m_nt->end();
      iter++){
    if((iter->getTag() & mask) == (tag & mask)){
      ((*iter).*mf)(arg);
    }
  }
}

/*!
  call a particle member function taking no argument for all particles
*/
template<typename T>
void ParallelParticleArray<T>::forAllParticles(void (T::*fnc)())
{
  for(typename NeighborTable<T>::iterator iter=m_nt->begin();
      iter!=m_nt->end();
      iter++){
    ((*iter).*fnc)();
  }
}

/*!
  call a const particle member function taking no argument for all particles
*/
template<typename T>
void ParallelParticleArray<T>::forAllParticles(void (T::*fnc)()const)
{
  for(typename NeighborTable<T>::iterator iter=m_nt->begin();
      iter!=m_nt->end();
      iter++){
    ((*iter).*fnc)();
  }
}

/*!
  call a particle member function taking one argument for all particles

  \param fnc the particle member function
  \param arg the argument to the particle member function
*/
template<typename T>
template<typename P>
void ParallelParticleArray<T>::forAllParticles(void (T::*fnc)(P),const P& arg)
{
  for(typename NeighborTable<T>::iterator iter=m_nt->begin();
      iter!=m_nt->end();
      iter++){
    ((*iter).*fnc)(arg);
  }
}

/*!
  call a particle member function taking one argument for all inner particles

  \param fnc the particle member function
  \param arg the argument to the particle member function
*/
template<typename T>
template<typename P>
void ParallelParticleArray<T>::forAllInnerParticles(void (T::*fnc)(P&),P& arg)
{

  NTBlock<T> InnerBlock=m_nt->inner();
  for(typename NTBlock<T>::iterator iter=InnerBlock.begin();
      iter!=InnerBlock.end();
      iter++){
    ((*iter).*fnc)(arg);
  }
}

/*!
  Call a constant particle member function taking no argument and returning a
  value for all particles and collect the return values in a container. The
  container has to be an STL sequence container (vector,list...) or something
  with the same interface. The template parameter P is a type of container of
  the return type of the particle member function, not the return type itself.
  The container had to be reference argument because template instantiation
  based only on return type is impossible.

  \param cont the container
  \param rdf the particle member function
*/
template<typename T>
template<typename P>
void ParallelParticleArray<T>::forAllParticlesGet(P& cont,typename P::value_type (T::*rdf)()const)
{
  for(typename NeighborTable<T>::iterator iter=m_nt->begin();
      iter!=m_nt->end();
      iter++){
    cont.push_back(((*iter).*rdf)());
  }
}

template<typename T>
ParallelParticleArray<T>::ParticleIterator::ParticleIterator(
  const NtBlock &ntBlock
)
  : m_ntBlock(ntBlock),
    m_it(m_ntBlock.begin())
{
  m_it = m_ntBlock.begin();
  m_numRemaining = m_ntBlock.size();
}

template<typename T>
bool ParallelParticleArray<T>::ParticleIterator::hasNext() const
{
  return (m_numRemaining > 0);
}

template<typename T>
typename ParallelParticleArray<T>::ParticleIterator::Particle &ParallelParticleArray<T>::ParticleIterator::next()
{
  Particle &p = (*m_it);
  m_it++;
  m_numRemaining--;
  return p;
}

template<typename T>
int ParallelParticleArray<T>::ParticleIterator::getNumRemaining() const
{
  return m_numRemaining;
}

template<typename T>
typename ParallelParticleArray<T>::ParticleIterator ParallelParticleArray<T>::getInnerParticleIterator()
{
  return ParticleIterator(m_nt->inner());
}

/*!
  Get a value for all inner particle using a particle member function and
  return the values in a container.

  \param cont the container
  \param rdf the particle member function
*/
template<typename T>
template<typename P>
void ParallelParticleArray<T>::forAllInnerParticlesGet(P& cont,typename P::value_type (T::*rdf)()const)
{
  NTBlock<T> InnerBlock=m_nt->inner();
  for(typename NTBlock<T>::iterator iter=InnerBlock.begin();
      iter!=InnerBlock.end();
      iter++){
    cont.push_back(((*iter).*rdf)());
  }
}

/*!
  Get a value for each particle using a particle member function and return a vector
  of pairs of the particle id and the value.

  \param rdf the particle member function
*/
template<typename T>
template <typename P>
vector<pair<int,P> > ParallelParticleArray<T>::forAllParticlesGetIndexed(P (T::*rdf)() const)
{
  vector<pair<int,P> > res;

  for(typename NeighborTable<T>::iterator iter=m_nt->begin();
      iter!=m_nt->end();
      iter++){
    res.push_back(make_pair(iter->getID(),((*iter).*rdf)()));
  }

  return res;
}

/*!
  Get a value all inner particles using a particle member function and return a vector
  of pairs of the particle id and the value.

  \param rdf the particle member function
*/
template<typename T>
template <typename P>
vector<pair<int,P> > ParallelParticleArray<T>::forAllInnerParticlesGetIndexed(P (T::*rdf)() const)
{
  vector<pair<int,P> > res;

  NTBlock<T> InnerBlock=m_nt->inner();
  int n=0;
  for(typename NTBlock<T>::iterator iter=InnerBlock.begin();
      iter!=InnerBlock.end();
      iter++){
    res.push_back(make_pair(iter->getID(),((*iter).*rdf)()));
    n++;
  }
  //cout<<" Inner particle number="<<n<<endl;
  return res;
}

/*!
  Call a constant particle member function taking no argument and returning a
  value for all particles which have a tag fitting a given tag and mask and
  collect the return values in a container.
  The container has to be an STL sequence container (vector,list...) or something
  with the same interface. The template parameter P is a type of container of
  the return type of the particle member function, not the return type itself.
  The container had to be reference argument because template instantiation
  based only on return type is impossible.

  \param cont the container
  \param rdf the particle member function
  \param tag the particle tag
  \param mask the mask
*/
template<typename T>
template<typename P>
void ParallelParticleArray<T>::forAllTaggedParticlesGet(P& cont,typename P::value_type (T::*rdf)()const,int tag,int mask)
{
  for(typename NeighborTable<T>::iterator iter=m_nt->begin();
      iter!=m_nt->end();
      iter++){
    if((iter->getTag() | mask )==(tag | mask)){
      cont.push_back(((*iter).*rdf)());
    }
  }
}

/*!
  Get a value for all inner particle which have a tag fitting a given tag and mask
  using a particle member function and return the values in a container.

  \param cont the container
  \param rdf the particle member function
  \param tag the particle tag
  \param mask the mask
*/
template<typename T>
template<typename P>
void ParallelParticleArray<T>::forAllTaggedInnerParticlesGet(P& cont,typename P::value_type (T::*rdf)()const,int tag,int mask)
{
  NTBlock<T> InnerBlock=m_nt->inner();
  for(typename NTBlock<T>::iterator iter=InnerBlock.begin();
      iter!=InnerBlock.end();
      iter++){
    int ptag=iter->getTag();
    if((ptag & mask )==(tag & mask)){
      cont.push_back(((*iter).*rdf)());
    }
  }
}

/*!
  Get a value for each particle which has a tag fitting a given tag and mask using
  a particle member function and return a vector of pairs of the particle id and the value.

  \param rdf the particle member function
  \param tag the particle tag
  \param mask the mask
*/
template<typename T>
template <typename P>
vector<pair<int,P> > ParallelParticleArray<T>::forAllTaggedParticlesGetIndexed(P (T::*rdf)() const,int tag,int mask)
{
  vector<pair<int,P> > res;

  for(typename NeighborTable<T>::iterator iter=m_nt->begin();
      iter!=m_nt->end();
      iter++){
    int id=iter->getID();
    int ptag=iter->getTag();
    if(( ptag & mask )==(tag & mask)){
      res.push_back(make_pair(id,((*iter).*rdf)()));
    }
  }

  return res;
}

/*!
  Get a value all inner particles which have a tag fitting a given tag and mask
  using a particle member function and return a vector of pairs of the particle
  id and the value.

  \param rdf the particle member function
  \param tag the particle tag
  \param mask the mask
*/
template<typename T>
template <typename P>
vector<pair<int,P> > ParallelParticleArray<T>::forAllInnerTaggedParticlesGetIndexed(P (T::*rdf)() const,int tag,int mask)
{
  vector<pair<int,P> > res;

  NTBlock<T> InnerBlock=m_nt->inner();
  for(typename NTBlock<T>::iterator iter=InnerBlock.begin();
      iter!=InnerBlock.end();
      iter++){
    int id=iter->getID();
    int ptag=iter->getTag();
    if((ptag & mask )==(tag & mask)){
      res.push_back(make_pair(id,((*iter).*rdf)()));
    }
  }

  return res;
}

/*!
  get a value for the nearest particle to each point of a grid using a particle member
  function and return a container with the values

  \param cont the container
  \param rdf the particle member function returning the value
  \param orig the base point of the grid
  \param dx the grid spacing in x-direction
  \param dy the grid spacing in y-direction
  \param dz the grid spacing in z-direction
  \param nx the grid size in x-direction
  \param ny the grid size in y-direction
  \param nz the grid size in z-direction
*/
template<typename T>
template <typename P>
void ParallelParticleArray<T>::forPointsGetNearest(P& cont,typename P::value_type (T::*rdf)() const,const Vec3& orig,
						   double dx,double dy,double dz,int nx,int ny,int nz)
{
  console.Debug() << "forPointsGetNearest" << "\n";

  Vec3 u_x=Vec3(1.0,0.0,0.0);
  Vec3 u_y=Vec3(0.0,1.0,0.0);
  Vec3 u_z=Vec3(0.0,0.0,1.0);
  for(int ix=0;ix<nx;ix++){
    for(int iy=0;iy<ny;iy++){
      for(int iz=0;iz<nz;iz++){
	Vec3 p=orig+double(ix)*dx*u_x+double(iy)*dy*u_y+double(iz)*dz*u_z;
	cont.push_back(((*(m_nt->getNearestPtr(p))).*rdf)());
      }
    }
  }

  console.Debug() << "end forPointsGetNearest" << "\n";
}

/*!
  Get the Ids of all particles in the boundary slab.

  \param dir the direction ,i.e. 0->x, 1->y and 2->z
  \param up up (1) or down (-1)
*/
template<typename T>
set<int> ParallelParticleArray<T>::getBoundarySlabIds(int dir,int up) const
{
  set<int> res;
  NTSlab<T> slab,slab2;

  switch(dir){
  case 0 :  // x-dir
    if(up==1){
      slab=m_nt->yz_slab(m_nt->xsize()-1);
      slab2=m_nt->yz_slab(m_nt->xsize()-2);
    } else {
      slab=m_nt->yz_slab(0);
      slab2=m_nt->yz_slab(1);
    }
    break;
  case 1 : // y-dir
    if(up==1){
      slab=m_nt->xz_slab(m_nt->ysize()-1);
      slab2=m_nt->xz_slab(m_nt->ysize()-2);
    } else {
      slab=m_nt->xz_slab(0);
      slab2=m_nt->xz_slab(1);
    }
    break;
  case 2 : // z-dir
    if(up==1){
      slab=m_nt->xy_slab(m_nt->zsize()-1);
      slab2=m_nt->xy_slab(m_nt->zsize()-2);
    } else {
      slab=m_nt->xy_slab(0);
      slab2=m_nt->xy_slab(1);
    }
    break;
  default:
    cout << "Error: wrong direction " << dir << " in getBoundarySlabIds" << endl;
  }

  for(typename NTSlab<T>::iterator iter=slab.begin();
      iter!=slab.end();
      iter++){
    res.insert(iter->getID());
  }
  for(typename NTSlab<T>::iterator iter=slab2.begin();
      iter!=slab2.end();
      iter++){
    res.insert(iter->getID());
  }

  return res;
}

/*!
  Get the Ids of all particles in the slab next to the boundary.

  \param dir the direction ,i.e. 0->x, 1->y and 2->z
  \param up up (1) or down (-1)
*/
template<typename T>
set<int> ParallelParticleArray<T>::get2ndSlabIds(int dir,int up) const
{
  set<int> res;
  NTSlab<T> slab;

  switch(dir){
  case 0 :  // x-dir
    if(up==1){
      slab=m_nt->yz_slab(m_nt->xsize()-2);
    } else {
      slab=m_nt->yz_slab(1);
    }
    break;
  case 1 : // y-dir
    if(up==1){
      slab=m_nt->xz_slab(m_nt->ysize()-2);
    } else {
      slab=m_nt->xz_slab(1);
    }
    break;
  case 2 : // z-dir
    if(up==1){
      slab=m_nt->xy_slab(m_nt->zsize()-2);
    } else {
      slab=m_nt->xy_slab(1);
    }
    break;
  default:
    cout << "Error: wrong direction " << dir << " in get2ndSlabIds" << endl;
  }

  for(typename NTSlab<T>::iterator iter=slab.begin();
      iter!=slab.end();
      iter++){
    res.insert(iter->getID());
  }

  return res;
}

/*!
  get all particles in inner block and put them into a vector

  \param pv a reference to the vector
*/
template<typename T>
void ParallelParticleArray<T>::getAllInnerParticles(vector<T>& pv)
{
  NTBlock<T> InnerBlock=m_nt->inner();
  for(typename NTBlock<T>::iterator iter=InnerBlock.begin();
      iter!=InnerBlock.end();
      iter++){
    pv.push_back(*iter);
  }
}


/*!
  save checkpoint data into an ostream

  \param ost the output stream
*/
template<typename T>
void ParallelParticleArray<T>::saveCheckPointData(std::ostream& ost)
{
  console.Debug() << "ParallelParticleArray<T>::saveCheckPointData\n";

  // output dimensions
  ost << m_nt->base_point() << "\n";
  ost << m_nt->base_idx_x() << " " << m_nt->base_idx_y() << " " << m_nt->base_idx_z() << "\n";

  // get nr. of particles in inner block
  ost << getInnerSize() << "\n";

  // save particles to stream
  NTBlock<T> InnerBlock=m_nt->inner();
  for(typename NTBlock<T>::iterator iter=InnerBlock.begin();
      iter!=InnerBlock.end();
      iter++){
    iter->saveCheckPointData(ost);
    ost << "\n";
  }
}

/*!
  load checkpoint data from an istream

  \param ist the input stream
*/
template<typename T>
void ParallelParticleArray<T>:: loadCheckPointData(std::istream& ist)
{
  console.Debug() << "ParallelParticleArray<T>::loadCheckPointData\n";
  int nparts;

  // get dimensions
  Vec3 bp;
  int bix,biy,biz;
  ist >> bp;
  ist >> bix;
  ist >> biy;
  ist >> biz;
  // barf if different from current values
  if((bp!=m_nt->base_point()) ||
     (bix!=m_nt->base_idx_x()) ||
     (biy!=m_nt->base_idx_y()) ||
     (biz!=m_nt->base_idx_z())){
    std::cerr << "local area data inconsistet: bp " << bp << " / " << m_nt->base_point()
	      << " bix: " << bix << " / " << m_nt->base_idx_x()
	      << " biy: " << biy << " / " << m_nt->base_idx_y()
	      << " bix: " << biz << " / " << m_nt->base_idx_z() << std::endl;
  }

  // get nr. of particles
  ist >> nparts;

  // load particles from stream and try to insert them
  for(int i=0;i!=nparts;i++){
    T p;
    p.loadCheckPointData(ist);
    if(!isInInner(p.getPos())) std::cerr << "not in inner: " << p.getPos() << std::endl;
    m_nt->insert(p);
  }
}


/*!
  output a particle array to an iostream. -- for debugging only --

  \param ost the output stream
  \param ppa the particle array
*/
template<typename T>
ostream& operator<<(ostream& ost, const ParallelParticleArray<T> &ppa)
{
  ost << "--ParallelParticleArray--" << endl;
  ost << *(ppa.m_nt) << endl;
  return ost;
}

