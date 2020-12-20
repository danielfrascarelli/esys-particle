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

#include "nt_block.h"

/*!
  construct an empty, uninitialized NeighborTable -> not usable
*/
template<typename T>
NeighborTable<T>::NeighborTable()
  : m_list(),
    m_array(),
    m_idParticleMap(),
    m_p0_global(),
    m_pmax_global(), //fluid contents
    m_dim(0),
    m_alpha(0),
    m_global_idx(0),
    m_global_idy(0),
    m_global_idz(0),
    m_xsize(0),
    m_ysize(0),
    m_zsize(0),
    m_valid(false),
    m_fluidexist(false) //fluid contents
{
}

/*!
  construct neigbortable with known dimensions

  \param x nr. of grid points in x-direction
  \param y nr. of grid points in y-direction
  \param z nr. of grid points in z-direction
  \param range grid spacing
  \param alpha pair search cutoff
  \param p0_global minimal corner (origin) of the global search space
  \param ix x-index of the local origin
  \param iy y-index of the local origin
  \param iz z-index of the local origin
*/
template<typename T>
NeighborTable<T>::NeighborTable(
  int x,
  int y,
  int z,
  double range,
  double alpha,
  const Vec3& p0_global,
  const Vec3& pmax_global,  //fluid contents
  int ix,
  int iy,
  int iz
)
  : m_p0_global(p0_global),
    m_pmax_global(pmax_global), //fluid contents
    m_dim(range),
    m_alpha(alpha),
    m_xsize(x),
    m_ysize(y),
    m_zsize(z),
    m_valid(true),
    m_fluidexist(false) //fluid contents
{
  m_global_idx=ix;
  m_global_idy=iy;
  m_global_idz=iz;
  m_min_corner=Vec3(m_p0_global.X()+m_dim*double(m_global_idx),
		    m_p0_global.Y()+m_dim*double(m_global_idy),
		    m_p0_global.Z()+m_dim*double(m_global_idz));
  m_max_corner=Vec3(m_p0_global.X()+m_dim*double(m_global_idx+m_xsize),
		    m_p0_global.Y()+m_dim*double(m_global_idy+m_ysize),
		    m_p0_global.Z()+m_dim*double(m_global_idz+m_zsize));
  console.Debug() << "Ntable corners : " << m_min_corner << " - " << m_max_corner << "\n";
  m_array.resize(m_xsize*m_ysize*m_zsize);
}

/*!
  destruct NeighborTable
*/
template<typename T>
NeighborTable<T>::~NeighborTable()
{
}

/*!
  clean up the search array
*/
template<typename T>
void NeighborTable<T>::clear_search_array()
{
  for(unsigned int iter=0;iter<m_array.size();iter++){
    m_array[iter].erase(m_array[iter].begin(),m_array[iter].end());
  }
}


//fluid contents
/*! 
  clean up the cell array
*/
template<typename T>
void NeighborTable<T>::clear_cell_array()
{
  for(unsigned int iter=0;iter<m_cellarray.size();iter++){
    m_cellarray[iter].erase(m_cellarray[iter].begin(),m_cellarray[iter].end());
  }
}


/*!
  helper function which returns the index in the search array
  into which a given position is mapped. Returns -1 if pos is
  outside the search space

  \param the position
*/
template<typename T>
int NeighborTable<T>::index(const Vec3& pos)
{
  int idx=-1;
  int ix=int(floor((pos.X()-m_p0_global.X())/m_dim))-m_global_idx;
  int iy=int(floor((pos.Y()-m_p0_global.Y())/m_dim))-m_global_idy;
  int iz=int(floor((pos.Z()-m_p0_global.Z())/m_dim))-m_global_idz;
  if((ix>=0)&&(ix<m_xsize)&&
     (iy>=0)&&(iy<m_ysize)&&
     (iz>=0)&&(iz<m_zsize)){
    idx=m_ysize*m_zsize*ix+m_zsize*iy+iz;
  }
  return idx;
}


//fluid contents
/*! 
  helper function which returns the index in the cell array
  into which a given position is mapped. Returns -1 if pos is
  outside the search space

  \param the position
*/
template<typename T>
int NeighborTable<T>::cellindex(const Vec3& pos)
{
  int idx=-1;

  int ix=int(floor((pos.X()-m_min_corner.X()-m_dim)/m_xside))+1;
  int iy=int(floor((pos.Y()-m_min_corner.Y()-m_dim)/m_yside))+1;
  int iz=int(floor((pos.Z()-m_min_corner.Z()-m_dim)/m_zside))+1;
  if((ix>=1)&&(ix<m_xcell-1)&&
     (iy>=1)&&(iy<m_ycell-1)&&
     (iz>=1)&&(iz<m_zcell-1)){
    idx=m_ycell*m_zcell*ix+m_zcell*iy+iz;
  }
  return idx;
}


/*!
  helper function which returns the index in the search array
  from given x-,y- and z-indices.

  \param x the x-index
  \param y the y-index
  \param z the z-index
  \warning no checks
*/
template<typename T>
int NeighborTable<T>::index(int x,int y,int z) const
{
  return m_ysize*m_zsize*x+m_zsize*y+z;
}


//fluid contents
/*!
  helper function which returns the index in the cell array
  from given x-,y- and z-indices.

  \param x the x-index
  \param y the y-index
  \param z the z-index
  \warning no checks
*/
template<typename T>
int NeighborTable<T>::cellindex(int x,int y,int z) const
{
  return m_ycell*m_zcell*x+m_zcell*y+z;
}


/*!
  check if a position is in the inner part

  \param pos the position
*/
template<typename T>
bool NeighborTable<T>::isInInner(const Vec3& pos)
{
  bool res;

  int ix=int(floor((pos.X()-m_p0_global.X())/m_dim))-m_global_idx;
  int iy=int(floor((pos.Y()-m_p0_global.Y())/m_dim))-m_global_idy;
  int iz=int(floor((pos.Z()-m_p0_global.Z())/m_dim))-m_global_idz;

  if(m_zsize<3) { // 2D-> ignore Z
    res=((ix>0) && (ix<m_xsize-1) && (iy>0) && (iy<m_ysize-1));
  } else {
    res=((ix>0) && (ix<m_xsize-1) && (iy>0) && (iy<m_ysize-1)&& (iz>0) && (iz<m_zsize-1));
  }

  return res;
}


/*!
  insert a particle into the NeighborTable

  \param t the particle
*/
template<typename T>
void NeighborTable<T>::insert(const T& t)
{
  // put into list
  typename list<T>::iterator iter=m_list.insert(m_list.end(),t);
  // put into array
  int idx=index(iter->getPos());
  if(idx!=-1){
//     m_idParticleMap.insert(IdParticleMap::value_type(iter->getID(),&(*iter)));
    m_idParticleMap.insert(make_pair(iter->getID(),&(*iter)));
    m_array[idx].push_back(iter);

    //fluid contents:
    if(m_fluidexist){
      int cellidx=cellindex(iter->getPos());
      if(cellidx!=-1){
        m_cellarray[cellidx].push_back(iter);
      };
    };

  } else { // outside -> delete
    typename list<T>::iterator h=iter;
    iter--;
    m_list.erase(h);
  }
}

/*!
  Build or rebuild the search array, calls clean_search_array.
  Particles outside the search space are removed from the list
*/
template<typename T>
void NeighborTable<T>::build()
{
  clear_search_array();
  if(m_fluidexist){clear_cell_array();}; //fluid contents
  // clean out id map
  m_idParticleMap.clear();
  for(typename list<T>::iterator iter=m_list.begin();
      iter!=m_list.end();
      iter++){
    int idx=index(iter->getPos());
    int id=iter->getID();
    if(idx!=-1){
      m_array[idx].push_back(iter);
      m_idParticleMap[id]=&(*iter);
      //fluid contents
      if(m_fluidexist){
        int cellidx=cellindex(iter->getPos());
        if(cellidx!=-1){m_cellarray[cellidx].push_back(iter);};
      }
    } else { // outside -> delete
      typename list<T>::iterator h=iter;
      iter--;
      // delete from list
      m_list.erase(h);
    }
  }
}

/*!
   Return pointer to particle at index.
*/
template<typename T>
T* NeighborTable<T>::ptr(NeighborTable<T>::indextype idx)
{
  return &(*(m_array[idx.first][idx.second]));
}

/*!
  Return reference to particle at index.
*/
template<typename T>
T& NeighborTable<T>::ref(NeighborTable<T>::indextype idx)
{
  return *(m_array[idx.first][idx.second]);
}

/*!
   Return pointer to particle with given id. Return NULL if the
   table doesn't contain the particle.

   \param id the id of the particle
*/
template<typename T>
T* NeighborTable<T>::ptr_by_id(int id)
{
  T* pParticle = NULL;
  typename IdParticleMap::iterator it = m_idParticleMap.find(id);

  if (it != m_idParticleMap.end()) {
    pParticle = it->second;
  }
  // debug output
  if(pParticle!=NULL){
    if(id!=pParticle->getID()){
      console.Debug() << "inconsistent idParticleMap: " << id << " vs. " << pParticle->getID() << "\n";
    }
  }
  return pParticle;
}

/*!
  Return pointer to particle nearest to given position. Returns
  NULL if position is outside search area.

  \param pos position
*/
template<typename T>
T *NeighborTable<T>::getNearestPtr(const Vec3& pos)
{
  T* res=NULL;

  // get grid index for pos
  int idx=index(pos);

  double dist=3.0*(m_dim*m_dim); // squared max. dist (grid diagonal)
  if(idx!=-1){
    for(typename pointtype::iterator iter=m_array[idx].begin();
	iter!=m_array[idx].end();
	iter++){
      double ndist=(pos-(*iter)->getPos()).norm2();
      res=(ndist<dist) ? &(**iter) : res;
      dist=(ndist<dist) ? ndist : dist;
    }
  }

  return res;
}

/*!
  erase particle (not meant to be used directly)
*/
template<typename T>
void NeighborTable<T>::erase(NeighborTable<T>::indextype idx)
{
  m_idParticleMap.erase(m_array[idx.first][idx.second]->getID());
  m_list.erase(m_array[idx.first][idx.second]);
  m_array[idx.first].erase(m_array[idx.first].begin()+idx.second);
}

/*!
  Return representation for a slab of the search array in the xy-plane.

  \param z position of the slab in z-direction
*/
template<typename T>
NTSlab<T> NeighborTable<T>::xy_slab(int z)
{
  return NTSlab<T>(this,DSlice(z,m_ysize,m_zsize,m_xsize,m_ysize*m_zsize));
}

/*!
  Return representation for a slab of the search array in the xz-plane.

  \param y position of the slab in y-direction
*/
template<typename T>
NTSlab<T> NeighborTable<T>::xz_slab(int y)
{
  return NTSlab<T>(this,DSlice(y*m_zsize,m_zsize,1,m_xsize,m_ysize*m_zsize));
}

/*!
  Return representation for a slab of the search array in the yz-plane.

  \param x position of the slab in x-direction
*/
template<typename T>
NTSlab<T> NeighborTable<T>::yz_slab(int x)
{
  return NTSlab<T>(this,DSlice(x*m_ysize*m_zsize,m_zsize,1,m_ysize,m_zsize));
}

/*
  Return block of the search array

  \param xmin minimum index in x-dimension
  \param xmax maximum index in x-dimension
  \param ymin minimum index in y-dimension
  \param ymax maximum index in y-dimension
  \param zmin minimum index in z-dimension
  \param zmax maximum index in z-dimension
*/
template<typename T>
NTBlock<T> NeighborTable<T>::block(int xmin,int xmax,int ymin,int ymax,int zmin,int zmax)
{
  return NTBlock<T>(this,xmin,xmax,ymin,ymax,zmin,zmax);
}

/*
  Return block of the search array.

  \param vmin minimum position
  \param vmax maximum position
*/
template<typename T>
NTBlock<T> NeighborTable<T>::block(const Vec3& vmin,const Vec3& vmax)
{
  // minimum corner index
  int xmin=int(floor((vmin.X()-m_p0_global.X())/m_dim))-m_global_idx;
  int ymin=int(floor((vmin.Y()-m_p0_global.Y())/m_dim))-m_global_idy;
  int zmin=int(floor((vmin.Z()-m_p0_global.Z())/m_dim))-m_global_idz;
  // check for minimum outside
  xmin=(xmin < 0) ? 0 : xmin;
  ymin=(ymin < 0) ? 0 : ymin;
  zmin=(zmin < 0) ? 0 : zmin;

  // maximum corner index
  int xmax=int(floor((vmax.X()-m_p0_global.X())/m_dim))-m_global_idx;
  int ymax=int(floor((vmax.Y()-m_p0_global.Y())/m_dim))-m_global_idy;
  int zmax=int(floor((vmax.Z()-m_p0_global.Z())/m_dim))-m_global_idz;
  // check for maximum outside
  xmax=(xmax < m_xsize) ? xmax : m_xsize-1;
  ymax=(ymax < m_ysize) ? ymax : m_ysize-1;
  zmax=(zmax < m_zsize) ? zmax : m_zsize-1;

  return NTBlock<T>(this,xmin,xmax,ymin,ymax,zmin,zmax);
}



/*!
  Return the inner region, i.e. excluding the boundary, as a NTBlock
*/
template<typename T>
NTBlock<T> NeighborTable<T>::inner()
{
  return NTBlock<T>(this,1,m_xsize-2,1,m_ysize-2,1,m_zsize-2);
}



/*!
  Add pairs containing the particles at two given gridpoints to list.

  \param list the pairlist
  \param idx1 the index of the first gridpoint
  \param idx2 the index of the second gridpoint
*/
template<typename T>
void NeighborTable<T>::addPairsToList(T_Handle<pairlist> list,int idx1,int idx2)
{

  for(typename pointtype::iterator iter=m_array[idx1].begin();
      iter!=m_array[idx1].end();
      iter++){
    for(typename pointtype::iterator iter2=m_array[idx2].begin();
	iter2!=m_array[idx2].end();
	iter2++){
      double dist2=((*iter)->getPos()-(*iter2)->getPos()).norm2();
      double dmax=(*iter)->getRad()+(*iter2)->getRad()+m_alpha;
      if(dist2<=(dmax*dmax)){
	if((*iter)->getID()<(*iter2)->getID()){
	  list->push_back(make_pair(&(**iter),&(**iter2)));
	} else {
	  list->push_back(make_pair(&(**iter2),&(**iter)));
	}
      }
    }
  }
}

/*!
  Add pairs containing the particles at one given gridpoint to list.

  \param list the pairlist
  \param idx the index of the  gridpoint
*/
template<typename T>
void NeighborTable<T>::addPairsToListLocal(T_Handle<pairlist> list,int idx)
{
  if(m_array[idx].size()>=2){ // at least 2 particles here, or no pairs
    for(typename pointtype::iterator iter=m_array[idx].begin();
	iter!=m_array[idx].end()-1;
	iter++){
      for(typename pointtype::iterator iter2=iter+1;
	  iter2!=m_array[idx].end();
	  iter2++){
	double dist2=((*iter)->getPos()-(*iter2)->getPos()).norm2();
	double dmax=(*iter)->getRad()+(*iter2)->getRad()+m_alpha;
	if(dist2<=(dmax*dmax)){
	  if((*iter)->getID()<(*iter2)->getID()){
	    list->push_back(make_pair(&(**iter),&(**iter2)));
	  } else {
	    list->push_back(make_pair(&(**iter2),&(**iter)));
	  }
	}
      }
    }
  }
}

/*!
  Add pairs containing at least one flagged particle at two
  given gridpoints to list.

  \param list the pairlist
  \param idx1 the index of the first gridpoint
  \param idx2 the index of the second gridpoint
*/
template<typename T>
void NeighborTable<T>::addPairsToListFlagged(T_Handle<pairlist> list,int idx1,int idx2)
{

  for(typename pointtype::iterator iter=m_array[idx1].begin();
      iter!=m_array[idx1].end();
      iter++){
    if((*iter)->isFlagged()){
      for(typename pointtype::iterator iter2=m_array[idx2].begin();
	  iter2!=m_array[idx2].end();
	  iter2++){
	if((*iter2)->isFlagged()){
	  if((*iter)->getID()<(*iter2)->getID()){
	    list->push_back(make_pair(&(**iter),&(**iter2)));
	  } else {
	    list->push_back(make_pair(&(**iter2),&(**iter)));
	  }
	}
      }
    }
  }
}

/*!
  Add pairs containing at least one flagged particle at one
  given gridpoint to list.

  \param list the pairlist
  \param idx the index of the  gridpoint
*/
template<typename T>
void NeighborTable<T>::addPairsToListLocalFlagged(T_Handle<pairlist> list,int idx)
{
  if(m_array[idx].size()>=2){ // at least 2 particles here, or no pairs
    for(typename pointtype::iterator iter=m_array[idx].begin();
	iter!=m_array[idx].end()-1;
	iter++){
      if((*iter)->isFlagged()){
	for(typename pointtype::iterator iter2=iter+1;
	    iter2!=m_array[idx].end();
	    iter2++){
	  if((*iter2)->isFlagged()){
	    if((*iter)->getID()<(*iter2)->getID()){
	      list->push_back(make_pair(&(**iter),&(**iter2)));
	    } else {
	      list->push_back(make_pair(&(**iter2),&(**iter)));
	    }
	  }
	}
      }
    }
  }
}


/*!
  Create a full list of pair of neighboring particles and return handle to it.
*/
template<typename T>
T_Handle<typename NeighborTable<T>::pairlist> NeighborTable<T>::getFullList()
{
  T_Handle<typename NeighborTable<T>::pairlist> list=new typename NeighborTable<T>::pairlist;

  for(int ix=0;ix<m_xsize;ix++){
    for(int iy=0;iy<m_ysize;iy++){
      for(int iz=0;iz<m_zsize;iz++){
	// add pairs within the gridpoint
	addPairsToListLocal(list,index(ix,iy,iz));
	// get search range, considering boundaries
	int xmax=(ix<m_xsize-1) ? ix+1 : ix;
	int ymax=(iy<m_ysize-1) ? iy+1 : iy;
	int zmax=(iz<m_zsize-1) ? iz+1 : iz;
	int xmin=(ix>0) ? ix-1 : ix;
	int ymin=(iy>0) ? iy-1 : iy;
	int zmin=(iz>0) ? iz-1 : iz;
	for(int i=xmin;i<=xmax;i++){
	  for(int j=ymin;j<=ymax;j++){
	    for(int k=zmin;k<=zmax;k++){
	      int idx1=index(ix,iy,iz);
	      int idx2=index(i,j,k);
	      if(idx2>idx1){
		addPairsToList(list,idx1,idx2);
	      }
	    }
	  }
	}
      }
    }
  }
  return list;
}

/*!
  Create a list of all pairs of neighboring particles involving a "flagged" particles
  and return handle to it.
*/
template<typename T>
T_Handle<typename NeighborTable<T>::pairlist> NeighborTable<T>::getNewList()
{
  T_Handle<typename NeighborTable<T>::pairlist> nlist=new typename NeighborTable<T>::pairlist;

  for(int ix=0;ix<m_xsize;ix++){
    for(int iy=0;iy<m_ysize;iy++){
      for(int iz=0;iz<m_zsize;iz++){
	// add pairs within the gridpoint
	addPairsToListLocalFlagged(nlist,index(ix,iy,iz));
	// get search range, considering boundaries
	int xmax=(ix<m_xsize-1) ? ix+1 : ix;
	int ymax=(iy<m_ysize-1) ? iy+1 : iy;
	int zmax=(iz<m_zsize-1) ? iz+1 : iz;
	int xmin=(ix>0) ? ix-1 : ix;
	int ymin=(iy>0) ? iy-1 : iy;
	int zmin=(iz>0) ? iz-1 : iz;
	for(int i=xmin;i<=xmax;i++){
	  for(int j=ymin;j<=ymax;j++){
	    for(int k=zmin;k<=zmax;k++){
	      int idx1=index(ix,iy,iz);
	      int idx2=index(i,j,k);
	      if(idx2>idx1){
		addPairsToListFlagged(nlist,idx1,idx2);
	      }
	    }
	  }
	}
      }
    }
  }
  for(typename list<T>::iterator iter=m_list.begin();
      iter!=m_list.end();
      iter++){
    iter->setFlag(false);
  }
  return nlist;
}

/*!
  Get list of all particles along a given plane. Naive implementation, i.e. check all particles
  for distance to plane. The plane is given by one point an the normal.

  \param orig The origin of the plane
  \param normal The normal of the plane
*/
template<typename T>
T_Handle<typename NeighborTable<T>::particlelist> NeighborTable<T>::getParticlesAtPlane(const Vec3& orig, const Vec3& normal)
{
  T_Handle<typename NeighborTable<T>::particlelist> nlist=new typename NeighborTable<T>::particlelist;

  console.Debug() << "NeighborTable<T>::getParticlesAtPlane: m_dim = " << m_dim << "\n";
  for(typename list<T>::iterator iter=m_list.begin();
      iter!=m_list.end();
      iter++){
    double dist=(iter->getPos()-orig)*normal;
    if(fabs(dist)<m_dim){
      nlist->push_back(&(*iter));
    }
  }
  console.Debug() << "NeighborTable<T>::getParticlesAtPlane: found = " << nlist->size() << "particles\n";

  return nlist;
}

/*!
  Get list of all particles near a given sphere.
  Naive implementation, i.e. check all particles for distance to sphere.

  \param centre The centre of the sphere
  \param radius The radius of the sphere
*/
template<typename T>
T_Handle<typename NeighborTable<T>::particlelist> NeighborTable<T>::getParticlesNearSphere(const Vec3& centre, const double& radius)
{
  T_Handle<typename NeighborTable<T>::particlelist> nlist=new typename NeighborTable<T>::particlelist;

  console.Debug() << "NeighborTable<T>::getParticlesNearSphere: m_dim = " << m_dim << "\n";
  for(typename list<T>::iterator iter=m_list.begin();
      iter!=m_list.end();
      iter++){
    double dist=(iter->getPos()-centre).norm() - radius;
    if(fabs(dist)<m_dim){
      nlist->push_back(&(*iter));
    }
  }
  console.Debug() << "NeighborTable<T>::getParticlesNearSphere: found = " << nlist->size() << "particles\n";

  return nlist;
}

/*!
  Get list of all particles near a given triangle. Checks all particles at grid within range
  of the bounding box of the triangle; (could be improved)

  \param T the triangle
*/
template<typename T>
T_Handle<typename NeighborTable<T>::particlelist> NeighborTable<T>::getParticlesNearTriangle(const Triangle& Tr)
{
  T_Handle<typename NeighborTable<T>::particlelist> nlist=new typename NeighborTable<T>::particlelist;

  // work out search area from bounding box of the triangle
  Vec3 v_min=Tr.getBoundingBoxMin()-Vec3(m_dim,m_dim,m_dim);
  Vec3 v_max=Tr.getBoundingBoxMax()+Vec3(m_dim,m_dim,m_dim);
  // check if completely outside
  if((v_min.X()<m_max_corner.X())&&(v_min.Y()<m_max_corner.Y())&&(v_min.Z()<m_max_corner.Z())&&
     (v_max.X()>m_min_corner.X())&&(v_max.Y()>m_min_corner.Y())&&(v_max.Z()>m_min_corner.Z())){

    NTBlock<T> TriangleBlock=block(v_min,v_max);
    for(typename NTBlock<T>::iterator iter=TriangleBlock.begin();
	iter!=TriangleBlock.end();
	iter++){
      if(Tr.sep(iter->getPos()) < iter->getRad()+m_alpha){
	nlist->push_back(&(*iter));
      }
    }
  }

  return nlist;
}

/*!
  Get list of all particles near a given edge. Checks all particles at grid within range
  of the bounding box of the triangle; (could be improved - search area via Bresenham)

  \param  E the Edge
*/
template<typename T>
T_Handle<typename NeighborTable<T>::particlelist> NeighborTable<T>::getParticlesNearEdge(const AEdge* E)
{
  T_Handle<typename NeighborTable<T>::particlelist> nlist=new typename NeighborTable<T>::particlelist;

  // get search area via bounding box
  Vec3 v_min=E->getBoundingBoxMin()-Vec3(m_dim,m_dim,m_dim);
  Vec3 v_max=E->getBoundingBoxMax()+Vec3(m_dim,m_dim,m_dim);
  // check if completely outside
  if((v_min.X()<m_max_corner.X())&&(v_min.Y()<m_max_corner.Y())&&(v_min.Z()<m_max_corner.Z())&&
     (v_max.X()>m_min_corner.X())&&(v_max.Y()>m_min_corner.Y())&&(v_max.Z()>m_min_corner.Z())){

    NTBlock<T> SearchBlock=block(v_min,v_max);

    for(typename NTBlock<T>::iterator iter=SearchBlock.begin();
	iter!=SearchBlock.end();
	iter++){
      if(E->sep(iter->getPos()) < iter->getRad()+m_alpha){
	nlist->push_back(&(*iter));
      }
    }
  }

  return nlist;
}

/*!
  Get list of all particles near a given point.

  \param  p the point
*/
template<typename T>
T_Handle<typename NeighborTable<T>::particlelist> NeighborTable<T>::getParticlesNearPoint(const Vec3& p)
{
  T_Handle<typename NeighborTable<T>::particlelist> nlist=new typename NeighborTable<T>::particlelist;

  // get search area via bounding box
  Vec3 v_min=p-Vec3(m_dim,m_dim,m_dim);
  Vec3 v_max=p+Vec3(m_dim,m_dim,m_dim);
  // check if completely outside
  if((v_min.X()<m_max_corner.X())&&(v_min.Y()<m_max_corner.Y())&&(v_min.Z()<m_max_corner.Z())&&
     (v_max.X()>m_min_corner.X())&&(v_max.Y()>m_min_corner.Y())&&(v_max.Z()>m_min_corner.Z())){

    NTBlock<T> SearchBlock=block(v_min,v_max);

    for(typename NTBlock<T>::iterator iter=SearchBlock.begin();
	iter!=SearchBlock.end();
	iter++){
      if((p-iter->getPos()).norm() < iter->getRad()+m_alpha){
	nlist->push_back(&(*iter));
      }
    }
  }

  return nlist;
}

/*!
  Get list of all particles
*/
template<typename T>
T_Handle<typename NeighborTable<T>::particlelist> NeighborTable<T>::getAllParticles()
{
  T_Handle<typename NeighborTable<T>::particlelist> nlist=new typename NeighborTable<T>::particlelist;

  for(typename list<T>::iterator iter=m_list.begin();
      iter!=m_list.end();
      iter++){
    nlist->push_back(&(*iter));
  }

  return nlist;
}

/*!
  output operator for NeighborTable

  \param ost the output stream
  \param NT the NeighborTable
*/
template<typename T>
std::ostream& operator<<(std::ostream& ost, const NeighborTable<T>& NT)
{
  ost << "---NeighborTable---" << endl;
  ost << "3d array dimensions (x,y,z),size : (" << NT.m_xsize << "," << NT.m_ysize<< "," << NT.m_zsize << "), " << NT.m_array.size() << endl;
  ost << "search range : " << NT.m_dim << endl;
  ost << "--list--" << endl;
  ost << (NT.m_list).size() << " elements" << endl;
  for(typename list<T>::const_iterator iter=(NT.m_list).begin();
      iter!=(NT.m_list).end();
      iter++){
    ost << iter->getID() << " " << iter->getPos() << endl;
  }
  ost << "---search array---" << endl;
  for(int ix=0;ix<NT.m_xsize;ix++){
    for(int iy=0;iy<NT.m_ysize;iy++){
      for(int iz=0;iz<NT.m_zsize;iz++){
	unsigned int np=NT.nparts_at_gridpoint(NT.index(ix,iy,iz));
	int idx=NT.index(ix,iy,iz);
	//ost << "---" << endl;
	ost << "(" << ix << "," << iy << "," << iz << ") , [" << idx << "], " << np << " : ";
	for(unsigned int i=0;i<np;i++){
	  ost << (NT.m_array[idx])[i]->getID() << " ";
	}
	ost << endl;
      }
    }
  }
  return ost;
}



/****fluid contents: begin****/

/*!
  Adding fluid
*/
template<typename T>
void NeighborTable<T>::initiateFluid(
  double Bw,
  double Bp,
  double Mu,
  double alpha,
  double flowrate,
  double pressure,
  Vec3 inflow,
  Vec3 outflow,
  int nx,
  int ny,
  int nz,
  int nx_min,
  int ny_min,
  int nz_min
)
{
  m_Bw=Bw;m_Bp=Bp;m_Mu=Mu;m_adjust=alpha;m_flowrate=flowrate;m_pressure=pressure;m_inflow=inflow;m_outflow=outflow;
  m_xcell=nx;m_ycell=ny;m_zcell=nz;
  m_global_cellidx=nx_min;m_global_cellidy=ny_min;m_global_cellidz=nz_min;
  m_xside=m_dim*double(m_xsize-2)/double(nx-2);
  m_yside=m_dim*double(m_ysize-2)/double(ny-2);
  m_zside=m_dim*double(m_zsize-2)/double(nz-2);
  /*
  cout << "neighbor table info: "<<endl;
  cout << "    Bw=" << m_Bw << " Bp=" <<m_Bp << " Mu=" << m_Mu << " m_adjust=" << m_adjust << endl;
  cout << "    m_flowrate=" << m_flowrate << " m_pressure=" <<m_pressure << " m_inflow=" << m_inflow << " m_outflow=" << m_outflow<<endl;
  cout << "    m_xcell=" << m_xcell <<" m_ycell="<<m_ycell<<" m_zcell="<<m_zcell<< endl;
  cout << "    m_xside=" << m_xside <<" m_yside="<<m_yside<<" m_zside="<<m_zside<< endl;
  */
  m_cellarray.resize(m_xcell*m_ycell*m_zcell);
  m_cells.resize(m_xcell*m_ycell*m_zcell);
  for(int i=0;i<m_xcell;i++){
    for(int j=0;j<m_ycell;j++){
      for(int k=0;k<m_zcell;k++){
        double global_x=m_min_corner.X()+m_dim+(i-0.5)*m_xside;
        double global_y=m_min_corner.Y()+m_dim+(j-0.5)*m_yside;
        double global_z=m_min_corner.Z()+m_dim+(k-0.5)*m_zside;
        Vec3 global_pos=Vec3(global_x,global_y,global_z);
        int idx=cellindex(i,j,k);
        m_cells[idx].setPos(global_pos);
        m_cells[idx].setSize(Vec3(m_xside,m_yside,m_zside));
        m_cells[idx].setIndex(Vec3(m_global_cellidx+i,m_global_cellidy+j,m_global_cellidz+k));
      }
    }
  }
  m_fluidexist=true;
}


/*!
  Create a full list of pair of particle and cell and return handle to it
*/
template<typename T>
T_Handle<typename NeighborTable<T>::particlecelllist> NeighborTable<T>::getParticleCellList()
{
  T_Handle<typename NeighborTable<T>::particlecelllist> nlist=new typename NeighborTable<T>::particlecelllist;

  int number=0;
  for(typename list<T>::iterator iter=m_list.begin();
      iter!=m_list.end();
      iter++){
    int idx=cellindex(iter->getPos());
    if(idx!=-1){
     nlist->push_back(make_pair(&(*iter),&m_cells[idx]));
     number++;
    }
  }
  return nlist;
}

/*!
  Create a new list of pair of particle and cell and return handle to it
*/
template<typename T>
T_Handle<typename NeighborTable<T>::particlecelllist> NeighborTable<T>::getNewParticleCellList()
{
  T_Handle<typename NeighborTable<T>::particlecelllist> nlist=new typename NeighborTable<T>::particlecelllist;
  for(typename list<T>::iterator iter=m_list.begin();
      iter!=m_list.end();
      iter++){
    if(iter->isFlagged()){
      int idx=cellindex(iter->getPos());
      if(idx!=-1){
        nlist->push_back(make_pair(&(*iter),&m_cells[idx]));
      }
    }
  }
  return nlist;
}


/*!
  Calculate initial fluid cell porosity using measurement sphere method
*/

template<typename T>
void NeighborTable<T>::calPorositySphere0()
{
  int number=0;
  int num_cells=0;
  double tot_d=0.0;
  double tot_r3=0.0;
  double r1=(m_xside+m_yside+m_zside)/6.0;
  double cellVolume=3.1415926*4.0/3.0*r1*r1*r1;
  for(int i=1;i<m_xcell-1;i++){
    for(int j=1;j<m_ycell-1;j++){
      for(int k=1;k<m_zcell-1;k++){
        int idx=cellindex(i,j,k);
        Vec3 pos=m_cells[idx].getPos();
        double volume=0;
        if(i==1||i==m_xcell-2||j==1||j==m_ycell-2||k==1||k==m_zcell-2){
          int ix=int(floor((pos.X()-m_p0_global.X())/m_dim))-m_global_idx;
          int iy=int(floor((pos.Y()-m_p0_global.Y())/m_dim))-m_global_idy;
          int iz=int(floor((pos.Z()-m_p0_global.Z())/m_dim))-m_global_idz;
          if((ix<1)||(ix>m_xsize-2)||(iy<1)||(iy>m_ysize-2)||(iz<1)||(iz>m_zsize-2)){
            cout<<"neighbor table size does not match fluid cell size, program aborted!"<<endl;
            abort();
          }
          double r=m_dim/2.0;
          double gridVolume=3.1415926*4.0/3.0*r*r*r;
          int ixmin=ix-1; int ixmax=ix+1;
          int iymin=iy-1; int iymax=iy+1;
          int izmin=iz-1; int izmax=iz+1;
          int inumber=0;
          double iPhi=0;
          for(int iix=ixmin;iix<=ixmax;iix++){
	    for(int iiy=iymin;iiy<=iymax;iiy++){
              for(int iiz=izmin;iiz<=izmax;iiz++){
                if((iix>=1)&&(iix<=m_xsize-2)&&(iiy>=1)&&(iiy<=m_ysize-2)&&(iiz>=1)&&(iiz<=m_zsize-2)){
                  double global_x=m_p0_global.X()+(iix+m_global_idx+0.5)*m_dim;
                  double global_y=m_p0_global.Y()+(iiy+m_global_idy+0.5)*m_dim;
                  double global_z=m_p0_global.Z()+(iiz+m_global_idz+0.5)*m_dim;
                  Vec3 nt_pos=Vec3(global_x,global_y,global_z);
                  int xmin=iix-1; int xmax=iix+1;
                  int ymin=iiy-1; int ymax=iiy+1;
                  int zmin=iiz-1; int zmax=iiz+1;
                  double ivolume=0;
                  for(int x=xmin;x<=xmax;x++){
	            for(int y=ymin;y<=ymax;y++){
	              for(int z=zmin;z<=zmax;z++){
                        int nt_idx=index(x,y,z);
                        for(typename pointtype::iterator iter=m_array[nt_idx].begin();
                            iter!=m_array[nt_idx].end();
                            iter++){
                          double r2=(*iter)->getRad();
                          double d=((*iter)->getPos()-nt_pos).norm();
                          if(d<=r-r2){
                            ivolume+=3.1415926*4.0/3.0*r2*r2*r2;
                          } else if(d<r+r2){
                            ivolume+=3.1415926/12.0/d*(r+r2-d)*(r+r2-d)*(d*d+2.0*d*(r+r2)-3.0*(r-r2)*(r-r2));
                          }
                        }
                      }
                    }
                  }
                  double phi=1.0-ivolume/gridVolume;
                  if(phi<=0){phi=0.00001;};
                  if(phi>=1){phi=0.99999;};
                  iPhi+=phi;
                  inumber++;
                }
              }
            }
          }

          double newPhi=iPhi/double(inumber);
          if(newPhi<=0.0){newPhi=0.00001;};
          if(newPhi>=1.0){newPhi=0.99999;};
          m_cells[idx].setPhi(newPhi);
          m_cells[idx].setnewPhi(newPhi);
          m_cells[idx].seteffPhi(newPhi);
        } else {
        int xmin=i-1; int xmax=i+1;
        int ymin=j-1; int ymax=j+1;
        int zmin=k-1; int zmax=k+1;
        for(int x=xmin;x<=xmax;x++){
	  for(int y=ymin;y<=ymax;y++){
	    for(int z=zmin;z<=zmax;z++){
	      int index=cellindex(x,y,z);
              for(typename pointtype::iterator iter=m_cellarray[index].begin();
                  iter!=m_cellarray[index].end();
                  iter++){
                double r2=(*iter)->getRad();
                double d=((*iter)->getPos()-pos).norm();
                if(d<=r1-r2){
                  volume+=3.1415926*4.0/3.0*r2*r2*r2;
                } else if(d<r1+r2){
                  volume+=3.1415926/12.0/d*(r1+r2-d)*(r1+r2-d)*(d*d+2.0*d*(r1+r2)-3.0*(r1-r2)*(r1-r2));
                }
              }
            }
          }
        }
        m_cells[idx].setVolume(volume);
        double newPhi=1.0-m_cells[idx].getVolume()/cellVolume;
        if(newPhi<=0.0){newPhi=0.00001;};
        if(newPhi>=1.0){newPhi=0.99999;};
        m_cells[idx].setPhi(newPhi);
        m_cells[idx].setnewPhi(newPhi);
        m_cells[idx].seteffPhi(newPhi);
        }
        int num_p=0;
        for(typename pointtype::iterator iter=m_cellarray[idx].begin();
            iter!=m_cellarray[idx].end();
            iter++) {
          double r=(*iter)->getRad();
          tot_d+=2.0*r;
          tot_r3+=pow(r,3.0);
          number++;
          num_p++;
        }
        if(num_p!=0){num_cells++;};
      }
    }
  }
  double meanK;
  if(number!=0){
    double meanD=tot_d/double(number);
    double tot_vp=4.0/3.0*3.1415926*tot_r3;
    double meanPhi=1.0-tot_vp/double(num_cells)/(m_xside*m_yside*m_zside);
    meanK=pow(meanD,2.0)/180.0*pow(meanPhi,3.0)/pow(1.0-meanPhi,2.0);
  } else{
    meanK=1e-7;
  }
  for(int i=1;i<m_xcell-1;i++){
    for(int j=1;j<m_ycell-1;j++){
      for(int k=1;k<m_zcell-1;k++){
        int idx=cellindex(i,j,k);
        double cell_totd=0.0;
        int p_in_cell=0;
        double newK;
        for(typename pointtype::iterator iter=m_cellarray[idx].begin();
            iter!=m_cellarray[idx].end();
            iter++) {
          double r=(*iter)->getRad();
          cell_totd+=2.0*r;
          p_in_cell++;
        }
        if(p_in_cell!=0){
          double D=cell_totd/double(p_in_cell);
          m_cells[idx].setD(D);//mean diameter of particles
          double Phi=m_cells[idx].getnewPhi();
          newK=pow(D,2.0)/180.0*pow(Phi,3.0)/pow(1.0-Phi,2.0);
          if(newK>meanK){newK=meanK;};
        } else {
          m_cells[idx].setD(0);
          newK=meanK;
        }
        m_cells[idx].setK(newK);
        m_cells[idx].seteffK(newK);
        double newBf=1.0/((1.0-m_cells[idx].getnewPhi())/m_Bp+m_cells[idx].getnewPhi()/m_Bw);
        m_cells[idx].setBf(newBf);
        m_cells[idx].seteffBf(newBf);
        m_cells[idx].setMu(m_Mu);
      }
    }
  }
}


/*!
  Calculate the fluid cell porosity using measurement sphere method
*/

template<typename T>
void NeighborTable<T>::calPorositySphere()
{
  double r1=(m_xside+m_yside+m_zside)/6.0;
  double cellVolume=3.1415926*4.0/3.0*r1*r1*r1;

  for(int i=1;i<m_xcell-1;i++){
    for(int j=1;j<m_ycell-1;j++){
      for(int k=1;k<m_zcell-1;k++){
        int idx=cellindex(i,j,k);
        Vec3 pos=m_cells[idx].getPos();
        double volume=0;

        if(i==1||i==m_xcell-2||j==1||j==m_ycell-2||k==1||k==m_zcell-2){
          int ix=int(floor((pos.X()-m_p0_global.X())/m_dim))-m_global_idx;
          int iy=int(floor((pos.Y()-m_p0_global.Y())/m_dim))-m_global_idy;
          int iz=int(floor((pos.Z()-m_p0_global.Z())/m_dim))-m_global_idz;
          if((ix<1)||(ix>m_xsize-2)||(iy<1)||(iy>m_ysize-2)||(iz<1)||(iz>m_zsize-2)){
            cerr<<"neighbor table size does not match fluid cell size, program aborted!"<<endl;
            abort();
          }
          double r=m_dim/2.0;
          double gridVolume=3.1415926*4.0/3.0*r*r*r;
          int ixmin=ix-1; int ixmax=ix+1;
          int iymin=iy-1; int iymax=iy+1;
          int izmin=iz-1; int izmax=iz+1;
          int inumber=0;
          double iPhi=0;
          for(int iix=ixmin;iix<=ixmax;iix++){
	        for(int iiy=iymin;iiy<=iymax;iiy++){
              for(int iiz=izmin;iiz<=izmax;iiz++){
                if((iix>=1)&&(iix<=m_xsize-2)&&(iiy>=1)&&(iiy<=m_ysize-2)&&(iiz>=1)&&(iiz<=m_zsize-2)){
                  double global_x=m_p0_global.X()+(iix+m_global_idx+0.5)*m_dim;
                  double global_y=m_p0_global.Y()+(iiy+m_global_idy+0.5)*m_dim;
                  double global_z=m_p0_global.Z()+(iiz+m_global_idz+0.5)*m_dim;
                  Vec3 nt_pos=Vec3(global_x,global_y,global_z);
                  int xmin=iix-1; int xmax=iix+1;
                  int ymin=iiy-1; int ymax=iiy+1;
                  int zmin=iiz-1; int zmax=iiz+1;
                  double ivolume=0;
                  for(int x=xmin;x<=xmax;x++){
	                for(int y=ymin;y<=ymax;y++){
	                  for(int z=zmin;z<=zmax;z++){
                        int nt_idx=index(x,y,z);
                        for(typename pointtype::iterator iter=m_array[nt_idx].begin();
                            iter!=m_array[nt_idx].end();
                            iter++){
                          double r2=(*iter)->getRad();
                          double d=((*iter)->getPos()-nt_pos).norm();
                          if(d<=r-r2){
                            ivolume+=3.1415926*4.0/3.0*r2*r2*r2;
                          } else if(d<r+r2){
                            ivolume+=3.1415926/12.0/d*(r+r2-d)*(r+r2-d)*(d*d+2.0*d*(r+r2)-3.0*(r-r2)*(r-r2));
                          }
                        }
                      }
                    }
                  }
                  double phi=1.0-ivolume/gridVolume;
                  if(phi<=0){phi=0.00001;};
                  if(phi>=1){phi=0.99999;};
                  iPhi+=phi;
                  inumber++;
                }
              }
            }
          }
          m_cells[idx].setPhi(m_cells[idx].getnewPhi());//porosity
          m_cells[idx].setnewPhi(iPhi/double(inumber));//new porosity
          if(m_cells[idx].getnewPhi()<=0){m_cells[idx].setnewPhi(0.00001);};
          if(m_cells[idx].getnewPhi()>=1){m_cells[idx].setnewPhi(0.99999);};
        } else {
        int xmin=i-1; int xmax=i+1;
        int ymin=j-1; int ymax=j+1;
        int zmin=k-1; int zmax=k+1;
        for(int x=xmin;x<=xmax;x++){
	  for(int y=ymin;y<=ymax;y++){
            for(int z=zmin;z<=zmax;z++){
	      int index=cellindex(x,y,z);
              for(typename pointtype::iterator iter=m_cellarray[index].begin();
                  iter!=m_cellarray[index].end();
                  iter++){
                double r2=(*iter)->getRad();
                double d=((*iter)->getPos()-pos).norm();
                if(d<=r1-r2){
                  volume+=3.1415926*4.0/3.0*r2*r2*r2;
                } else if(d<r1+r2){
                  volume+=3.1415926/12.0/d*(r1+r2-d)*(r1+r2-d)*(d*d+2.0*d*(r1+r2)-3.0*(r1-r2)*(r1-r2));
                }
              }
            }
          }
        }
        m_cells[idx].setVolume(volume);
        m_cells[idx].setPhi(m_cells[idx].getnewPhi());//porosity
        m_cells[idx].setnewPhi(1.0-m_cells[idx].getVolume()/cellVolume);//new porosity
        if(m_cells[idx].getnewPhi()<=0){m_cells[idx].setnewPhi(0.00001);};
        if(m_cells[idx].getnewPhi()>=1){m_cells[idx].setnewPhi(0.99999);};
        }

        m_cells[idx].seteffPhi((1.0-m_adjust)*m_cells[idx].getPhi()+m_adjust*m_cells[idx].getnewPhi());//effective porosity
      }
    }
  }
}


/*!
  Update fluid cell: m_D, m_Vp, m_K, m_Bf;
*/
template<typename T>
void NeighborTable<T>::updateFluidcells(){
  int number=0;
  int num_cells=0;
  double tot_d=0.0;
  double tot_r3=0.0;
  for(int i=1;i<m_xcell-1;i++){
    for(int j=1;j<m_ycell-1;j++){
      for(int k=1;k<m_zcell-1;k++){
        int idx=cellindex(i,j,k);
        int num_p=0;
        for(typename pointtype::iterator iter=m_cellarray[idx].begin();
            iter!=m_cellarray[idx].end();
            iter++) {
          double r=(*iter)->getRad();
          tot_d+=2.0*r;
          tot_r3+=pow(r,3.0);
          number++;
          num_p++;
        }
        if(num_p!=0){num_cells++;};        
      }
    }
  }

  double meanK;
  if(number!=0){
    double meanD=tot_d/double(number);
    double tot_vp=4.0/3.0*3.1415926*tot_r3;
    double meanPhi=1.0-tot_vp/double(num_cells)/(m_xside*m_yside*m_zside);
    meanK=pow(meanD,2.0)/180.0*pow(meanPhi,3.0)/pow(1.0-meanPhi,2.0);
  } else{
    meanK=1e-7;
  }  
  for(int i=1;i<m_xcell-1;i++){
    for(int j=1;j<m_ycell-1;j++){
      for(int k=1;k<m_zcell-1;k++){
        int idx=cellindex(i,j,k);
        double cell_totd=0.0;
        int p_in_cell=0;
        double newK;
        Vec3 tot_vel;
        for(typename pointtype::iterator iter=m_cellarray[idx].begin();
            iter!=m_cellarray[idx].end();
            iter++) {
          double r=(*iter)->getRad();
          cell_totd+=2.0*r;
          p_in_cell++;
        }
        if(p_in_cell!=0){
          double D=cell_totd/double(p_in_cell);
          m_cells[idx].setD(D);//mean diameter of particles
          Vec3 Vp=tot_vel/double(p_in_cell);
          m_cells[idx].setVp(Vp);//mean particle velocity
          double Phi=m_cells[idx].getnewPhi();
          newK=pow(D,2.0)/180.0*pow(Phi,3.0)/pow(1.0-Phi,2.0);
          if(newK>meanK){newK=meanK;};
        } else {
          m_cells[idx].setD(0);
          m_cells[idx].setVp(Vec3(0.0,0.0,0.0));
          newK=meanK;
        }
        m_cells[idx].seteffK((1.0-m_adjust)*m_cells[idx].getK()+m_adjust*newK); //effective permeability
        m_cells[idx].setK(newK);//updated permeability
        double Phi=m_cells[idx].getnewPhi();
        double newBf=1.0/((1.0-Phi)/m_Bp+Phi/m_Bw);
        m_cells[idx].seteffBf((1.0-m_adjust)*m_cells[idx].getBf()+m_adjust*newBf);//effective bulk mudulus
        m_cells[idx].setBf(newBf);//updated bulk mudulus
        m_cells[idx].setMu(m_Mu);
      }
    }
  }
}


/*!
  Calculate coefficients for linear equations
*/
template<typename T>
void NeighborTable<T>::calCoeffi(double timestep,int nt) {

  double A_w, A_e, A_n, A_s, A_u, A_d, A_c, B_w, B_e, B_n, B_s, B_u, B_d, B_c;
  double co_w,co_e,co_n,co_s,co_u, co_d, co_c,co_b;
  int idx,idx_w,idx_e,idx_s,idx_n,idx_u,idx_d;
  double K,Kw,Ke,Ks,Kn,Ku,Kd;
  double P,Pw,Pe,Ps,Pn,Pu,Pd;
  double flowrate;
  double V=m_xside*m_yside*m_zside;
  for(int ix=1;ix<m_xcell-1;ix++){
    for(int iy=1;iy<m_ycell-1;iy++){
      for(int iz=1;iz<m_zcell-1;iz++){
        idx=cellindex(ix,iy,iz);
        idx_w=cellindex(ix-1,iy,iz); idx_e=cellindex(ix+1,iy,iz); idx_n=cellindex(ix,iy-1,iz);
        idx_s=cellindex(ix,iy+1,iz); idx_d=cellindex(ix,iy,iz-1); idx_u=cellindex(ix,iy,iz+1);
        K=m_cells[idx].geteffK();
        Kw=m_cells[idx_w].geteffK(); Ke=m_cells[idx_e].geteffK(); Ks=m_cells[idx_s].geteffK();
        Kn=m_cells[idx_n].geteffK(); Ku=m_cells[idx_u].geteffK(); Kd=m_cells[idx_d].geteffK();
        P=m_cells[idx].getdisP();
        Pw=m_cells[idx_w].getdisP(); Pe=m_cells[idx_e].getdisP(); Ps=m_cells[idx_s].getdisP();
        Pn=m_cells[idx_n].getdisP(); Pu=m_cells[idx_u].getdisP(); Pd=m_cells[idx_d].getdisP();
        double Bf=m_cells[idx].geteffBf();
        double Phi=m_cells[idx].geteffPhi();
        double dV=(m_cells[idx].getPhi()-m_cells[idx].getnewPhi())/nt; //pore volume change
        double dP=Bf*dV; //pressure change
        if(fabs(m_min_corner.X()+m_dim-m_p0_global.X())<1e-6 && ix==1){ //cell located at west boundary
          if(m_inflow.X()==1 || m_inflow.X()==2){ //inlet
            A_w=0;
            co_w=0;
            if(m_flowrate==0){flowrate=K*m_pressure/m_Mu/m_xside;}
            else{flowrate=m_flowrate;};
            B_w=Bf/V/Phi*flowrate*timestep*m_yside*m_zside;
          } else if(m_outflow.X()==-1 || m_outflow.X()==2){ //outlet
            A_w=K*m_xside*timestep/m_Mu;
            co_w=0;
            B_w=0;
          } else { //closed boundary
            A_w=0;co_w=0;B_w=0;
          }
        } else{ //cell located in inner area
          A_w=0.5*(K+Kw)*timestep*m_yside*m_zside/m_xside/m_Mu;
          co_w=Bf/V/Phi*A_w;
          B_w=(1.0-m_adjust)*co_w*Pw;
        }

        if(fabs(m_max_corner.X()-m_dim-m_pmax_global.X())<1e-6 && ix==(m_xcell-2)){ //cell located at east boundary
          if(m_inflow.X()==-1 || m_inflow.X()==2){ //inlet
            A_e=0;
            co_e=0;
            if(m_flowrate==0){flowrate=K*m_pressure/m_Mu/m_xside;}
            else{flowrate=m_flowrate;};
            B_e=Bf/V/Phi*flowrate*timestep*m_yside*m_zside;
          } else if(m_outflow.X()==1 || m_outflow.X()==2){ //outlet
            A_e=K*m_xside*timestep/m_Mu;
            co_e=0;
            B_e=0;
          } else { //closed boundary
            A_e=0;co_e=0;B_e=0;
          }
        } else{ //cell located in inner area
          A_e=0.5*(K+Ke)*timestep*m_yside*m_zside/m_xside/m_Mu;
          co_e=Bf/V/Phi*A_e;
          B_e=(1.0-m_adjust)*co_e*Pe;
        }

        if(fabs(m_min_corner.Y()+m_dim-m_p0_global.Y())<1e-6 && iy==1){ //cell located at north boundary
          if(m_inflow.Y()==1 || m_inflow.Y()==2){ //inlet
            A_n=0;
            co_n=0;
            if(m_flowrate==0){flowrate=K*m_pressure/m_Mu/m_yside;}
            else{flowrate=m_flowrate;};
            B_n=Bf/V/Phi*flowrate*timestep*m_xside*m_zside;
          } else if(m_outflow.Y()==-1 || m_outflow.Y()==2){ //outlet
            A_n=K*m_yside*timestep/m_Mu;
            co_n=0;
            B_n=0;
          } else { //closed boundary
            A_n=0;co_n=0;B_n=0;
          }
        } else{ //cell located in inner area
          A_n=0.5*(K+Kn)*timestep*m_xside*m_zside/m_yside/m_Mu;
          co_n=Bf/V/Phi*A_n;
          B_n=(1.0-m_adjust)*co_n*Pn;
        }

        if(fabs(m_max_corner.Y()-m_dim-m_pmax_global.Y())<1e-6 && iy==(m_ycell-2)){ //cell located at south boundary
          if(m_inflow.Y()==-1 || m_inflow.Y()==2){ //inlet
            A_s=0;
            co_s=0;
            if(m_flowrate==0){flowrate=K*m_pressure/m_Mu/m_yside;}
            else{flowrate=m_flowrate;};
            B_s=Bf/V/Phi*flowrate*timestep*m_xside*m_zside;
          } else if(m_outflow.Y()==1 || m_outflow.Y()==2){ //outlet
            A_s=K*m_yside*timestep/m_Mu;
            co_s=0;
            B_s=0;
          } else { //closed boundary
            A_s=0;co_s=0;B_s=0;
          }
        } else{ //cell located in inner area
          A_s=0.5*(K+Ks)*timestep*m_xside*m_zside/m_yside/m_Mu;
          co_s=Bf/V/Phi*A_s;
          B_s=(1.0-m_adjust)*co_s*Ps;
        }

       if(fabs(m_min_corner.Z()+m_dim-m_p0_global.Z())<1e-6 && iz==1){ //cell located at down/lower boundary
          if(m_inflow.Z()==1 || m_inflow.Z()==2){ //inlet
            A_d=0;
            co_d=0;
            if(m_flowrate==0){flowrate=K*m_pressure/m_Mu/m_zside;}
            else{flowrate=m_flowrate;};
            B_d=Bf/V/Phi*flowrate*timestep*m_xside*m_yside;
          } else if(m_outflow.Z()==-1 || m_outflow.Z()==2){ //outlet
            A_d=K*m_zside*timestep/m_Mu;
            co_d=0;
            B_d=0;
          } else { //closed boundary
            A_d=0;co_d=0;B_d=0;
          }
        } else{ //cell located in inner area
          A_d=0.5*(K+Kd)*timestep*m_xside*m_yside/m_zside/m_Mu;
          co_d=Bf/V/Phi*A_d;
          B_d=(1.0-m_adjust)*co_d*Pd;
        }

        if(fabs(m_max_corner.Z()-m_dim-m_pmax_global.Z())<1e-6 && iz==(m_zcell-2)){ //cell located at upper boundary
          if(m_inflow.Z()==-1 || m_inflow.Z()==2){ //inlet
            A_u=0;
            co_u=0;
            if(m_flowrate==0){flowrate=K*m_pressure/m_Mu/m_zside;}
            else{flowrate=m_flowrate;};
            B_u=Bf/V/Phi*flowrate*timestep*m_xside*m_yside;
          } else if(m_outflow.Z()==1 || m_outflow.Z()==2){ //outlet
            A_u=K*m_zside*timestep/m_Mu;
            co_u=0;
            B_u=0;
          } else { //closed boundary
            A_u=0;co_u=0;B_u=0;
          }
        } else{ //cell located in inner area
          A_u=0.5*(K+Ku)*timestep*m_xside*m_yside/m_zside/m_Mu;
          co_u=Bf/V/Phi*A_u;
          B_u=(1.0-m_adjust)*co_u*Pu;
        }
        A_c=-(A_w+A_e+A_n+A_s+A_d+A_u);
        co_c=m_adjust*Bf/V/Phi*A_c-1.0;
        B_c=((1.0-m_adjust)*A_c*Bf/V/Phi+1.0)*P; //coefficent for central cell
        co_b=-(B_w+B_e+B_n+B_s+B_d+B_u+B_c+dP); //Constant
        m_cells[idx].setc_W(m_adjust*co_w);
        m_cells[idx].setc_E(m_adjust*co_e);
        m_cells[idx].setc_N(m_adjust*co_n);
        m_cells[idx].setc_S(m_adjust*co_s);
        m_cells[idx].setc_D(m_adjust*co_d);
        m_cells[idx].setc_U(m_adjust*co_u);
        m_cells[idx].setc_C(co_c);
        m_cells[idx].setc_B(co_b);
      }
    }
  }
}


/*!
  setting distributed pore pressures;
*/
template<typename T>
void NeighborTable<T>::setPressure(vector<pair<Vec3,double> > pressure)
{
  for(int ix=1;ix<m_xcell-1;ix++){
    for(int iy=1;iy<m_ycell-1;iy++){
      for(int iz=1;iz<m_zcell-1;iz++){
        int idx=cellindex(ix,iy,iz);
        Vec3 global_index=m_cells[idx].getIndex();
        for(vector<pair<Vec3,double> >::iterator iter=pressure.begin();
            iter!=pressure.end();
            iter++){
          if(global_index==iter->first){
            
            if(iter->second<0){
              m_cells[idx].setdisP(0);
            } else {
              m_cells[idx].setdisP(iter->second);
            };
            break;
          };
        }
      }
    }
  }
  calVelocity();
}


/*!
  calculate fluid velocity and pressure gradient
*/
template<typename T>
void NeighborTable<T>::calVelocity()
{
  int idx,idx_w,idx_e,idx_s,idx_n,idx_u,idx_d;
  double K,Kw,Ke,Ks,Kn,Ku,Kd;
  double P,Pw,Pe,Ps,Pn,Pu,Pd;
  double Vw,Ve,Vn,Vs,Vd,Vu;
  double flowrate;

  for(int ix=1;ix<m_xcell-1;ix++){
    for(int iy=1;iy<m_ycell-1;iy++){
      for(int iz=1;iz<m_zcell-1;iz++){
        idx=cellindex(ix,iy,iz);
        idx_w=cellindex(ix-1,iy,iz); idx_e=cellindex(ix+1,iy,iz); idx_n=cellindex(ix,iy-1,iz);
        idx_s=cellindex(ix,iy+1,iz); idx_d=cellindex(ix,iy,iz-1); idx_u=cellindex(ix,iy,iz+1);
        K=m_cells[idx].getK();
        Kw=m_cells[idx_w].getK(); Ke=m_cells[idx_e].getK(); Ks=m_cells[idx_s].getK();
        Kn=m_cells[idx_n].getK(); Ku=m_cells[idx_u].getK(); Kd=m_cells[idx_d].getK();
        P=m_cells[idx].getdisP();
        Pw=m_cells[idx_w].getdisP(); Pe=m_cells[idx_e].getdisP(); Ps=m_cells[idx_s].getdisP();
        Pn=m_cells[idx_n].getdisP(); Pu=m_cells[idx_u].getdisP(); Pd=m_cells[idx_d].getdisP();
        double Phi=m_cells[idx].getnewPhi();
        //calculate fluid velocity:
        if(fabs(m_min_corner.X()+m_dim-m_p0_global.X())<1e-6 && ix==1){ //cell located at west boundary
          if(m_inflow.X()==1 || m_inflow.X()==2){
            if(m_flowrate==0){flowrate=K*m_pressure/m_Mu/m_xside;}
            else{flowrate=m_flowrate;};
            Vw=flowrate/Phi;Pw=flowrate*m_Mu*m_xside/K+P;
          } //inlet
          else if(m_outflow.X()==-1 || m_outflow.X()==2){
            Vw=K*(0-P)/m_xside/m_Mu/Phi;Pw=0;
          } //outlet
          else {Vw=999999;Pw=999999;} //closed boundary
        } else {
          Vw=0.5*(K+Kw)*(P-Pw)/m_xside/m_Mu/Phi;
        }

        if(fabs(m_max_corner.X()-m_dim-m_pmax_global.X())<1e-6 && ix==(m_xcell-2)){ //cell located at east boundary
          if(m_inflow.X()==-1 || m_inflow.X()==2){
            if(m_flowrate==0){flowrate=K*m_pressure/m_Mu/m_xside;}
            else{flowrate=m_flowrate;};
            Ve=flowrate/Phi;Pe=flowrate*m_Mu*m_xside/K+P;
          } //inlet
          else if(m_outflow.X()==1 || m_outflow.X()==2){
            Ve=K*(0-P)/m_xside/m_Mu/Phi;Pe=0;
          } //outlet
          else {Ve=999999;Pe=999999;} //closed boundary
        } else {
          Ve=0.5*(K+Ke)*(P-Pe)/m_xside/m_Mu/Phi;
        }

        if(fabs(m_min_corner.Y()+m_dim-m_p0_global.Y())<1e-6 && iy==1){ //cell located at north boundary
          if(m_inflow.Y()==1 || m_inflow.Y()==2){
            if(m_flowrate==0){flowrate=K*m_pressure/m_Mu/m_yside;}
            else{flowrate=m_flowrate;};
            Vn=flowrate/Phi;Pn=flowrate*m_Mu*m_yside/K+P;
          } //inlet
          else if(m_outflow.Y()==-1 || m_outflow.Y()==2){
            Vn=K*(0-P)/m_yside/m_Mu/Phi;Pn=0;
          } //outlet
          else {Vn=999999;Pn=999999;} //closed boundary
        } else {
          Vn=0.5*(K+Kn)*(P-Pn)/m_yside/m_Mu/Phi;
        }

        if(fabs(m_max_corner.Y()-m_dim-m_pmax_global.Y())<1e-6 && iy==(m_ycell-2)){ //cell located at south boundary
          if(m_inflow.Y()==-1 || m_inflow.Y()==2){
            if(m_flowrate==0){flowrate=K*m_pressure/m_Mu/m_yside;}
            else{flowrate=m_flowrate;};
            Vs=flowrate/Phi;Ps=flowrate*m_Mu*m_yside/K+P;
          } //inlet
          else if(m_outflow.Y()==1 || m_outflow.Y()==2){
            Vs=K*(0-P)/m_yside/m_Mu/Phi;Ps=0;
          } //outlet
          else {Vs=999999;Ps=999999;} //closed boundary
        } else {
          Vs=0.5*(K+Ks)*(P-Ps)/m_yside/m_Mu/Phi;
        }

        if(fabs(m_min_corner.Z()+m_dim-m_p0_global.Z())<1e-6 && iz==1){ //cell located at lower boundary
          if(m_inflow.Z()==1 || m_inflow.Z()==2){
            if(m_flowrate==0){flowrate=K*m_pressure/m_Mu/m_zside;}
            else{flowrate=m_flowrate;};
            Vd=flowrate/Phi;Pd=flowrate*m_Mu*m_zside/K+P;
          } //inlet
          else if(m_outflow.Z()==-1 || m_outflow.Z()==2){
            Vd=K*(0-P)/m_zside/m_Mu/Phi;Pd=0;
          } //outlet
          else {Vd=999999;Pd=999999;} //closed boundary
        } else {
          Vd=0.5*(K+Kd)*(P-Pd)/m_zside/m_Mu/Phi;
        }

        if(fabs(m_max_corner.Z()-m_dim-m_pmax_global.Z())<1e-6 && iz==(m_zcell-2)){ //cell located at upper boundary
          if(m_inflow.Z()==-1 || m_inflow.Z()==2){
            if(m_flowrate==0){flowrate=K*m_pressure/m_Mu/m_zside;}
            else{flowrate=m_flowrate;};
            Vu=flowrate/Phi;Pu=flowrate*m_Mu*m_zside/K+P;
          } //inlet
          else if(m_outflow.Z()==1 || m_outflow.Z()==2){
            Vu=K*(0-P)/m_zside/m_Mu/Phi;Pu=0;
          } //outlet
          else {Vu=999999;Pu=999999;} //closed boundary
        } else {
          Vu=0.5*(K+Ku)*(P-Pu)/m_zside/m_Mu/Phi;
        }

        double Vx,Vy,Vz,Px,Py,Pz;
        if(Vw==999999 && Ve==999999){
          Vx=0;Px=0;
        } else if(Vw==999999){
          Vx=Ve;Px=-Pe/m_xside;
        } else if(Ve==999999){
          Vx=-Vw;Px=Pw/m_xside;
        } else {
          Vx=(-Vw+Ve)/2.0;Px=(Pw-Pe)/2.0/m_xside;
        }

        if(Vn==999999 && Vs==999999){
          Vy=0;Py=0;
        } else if(Vn==999999){
          Vy=Vs;Py=-Ps/m_yside;
        } else if(Vs==999999){
          Vy=-Vn;Py=Pn/m_yside;
        } else {
          Vy=(-Vn+Vs)/2.0;Py=(Pn-Ps)/2.0/m_yside;
        }

        if(Vd==999999 && Vu==999999){
          Vz=0;Pz=0;
        } else if(Vd==999999){
          Vz=Vu;Pz=-Pu/m_zside;
        } else if(Vu==999999){
          Vz=-Vd;Pz=Pd/m_zside;
        } else {
          Vz=(-Vd+Vu)/2.0;Pz=(Pd-Pu)/2.0/m_zside;
        }
        m_cells[idx].setVf(Vec3(Vx,Vy,Vz));//fluid velocity
        m_cells[idx].setPg(Vec3(Px,Py,Pz));//pressure gradient
      }
    }
  }
}


template<typename T>
vector<CFluidCell> NeighborTable<T>::get_yz_cellslab(int x)
{
  vector<CFluidCell> res;
  for(int iy=1;iy<m_ycell-1;iy++){
    for(int iz=1;iz<m_zcell-1;iz++){
      int idx=cellindex(x,iy,iz);
      res.push_back(m_cells[idx]);
    }
  }
  return res;
}


template<typename T>
vector<CFluidCell> NeighborTable<T>::get_xz_cellslab(int y)
{
  vector<CFluidCell> res;;
  for(int ix=1;ix<m_xcell-1;ix++){
    for(int iz=1;iz<m_zcell-1;iz++){
      int idx=cellindex(ix,y,iz);
      res.push_back(m_cells[idx]);
    }
  }
  return res;
}


template<typename T>
vector<CFluidCell> NeighborTable<T>::get_xy_cellslab(int z)
{
  vector<CFluidCell> res;
  for(int ix=1;ix<m_xcell-1;ix++){
    for(int iy=1;iy<m_ycell-1;iy++){
      int idx=cellindex(ix,iy,z);
      res.push_back(m_cells[idx]);
    }
  }
  return res;
}



template<typename T>
void NeighborTable<T>::set_yz_cellslab(int x, vector<CFluidCell> recv_slab)
{
  if(recv_slab.size()==(m_ycell-2)*(m_zcell-2)){
    int number=0;
    for(int iy=1;iy<m_ycell-1;iy++){
      for(int iz=1;iz<m_zcell-1;iz++){
        int idx=cellindex(x,iy,iz);
        m_cells[idx]=recv_slab[number];
        number++;
      }
    }
  } else{
   //cerr << "FluidcellTable<T>::set_yz_slab ERROR: size mismatch in received data, recv_slab.size()!=m_ycell*m_zcell - (" << recv_slab.size() << "!=" << m_ycell*m_zcell << ")"<< std::endl;
  }
}

template<typename T>
void NeighborTable<T>::set_xz_cellslab(int y, vector<CFluidCell> recv_slab)
{
  if(recv_slab.size()==(m_xcell-2)*(m_zcell-2)){
    int number=0;
    for(int ix=1;ix<m_ycell-1;ix++){
      for(int iz=1;iz<m_zcell-1;iz++){
        int idx=cellindex(ix,y,iz);
        m_cells[idx]=recv_slab[number];
        number++;
      }
    }
  } else{
    //cerr << "FluidcellTable<T>::set_xz_slab ERROR: size mismatch in received data, recv_slab.size()!=m_xcell*m_zcell - (" << recv_slab.size() << "!=" << m_xcell*m_zcell << ")"<< std::endl;
  }
}

template<typename T>
void NeighborTable<T>::set_xy_cellslab(int z, vector<CFluidCell> recv_slab)
{
  if(recv_slab.size()==(m_xcell-2)*(m_ycell-2)){
    int number=0;
    for(int ix=1;ix<m_xcell-1;ix++){
      for(int iy=1;iy<m_ycell-1;iy++){
        int idx=cellindex(ix,iy,z);
        m_cells[idx]=recv_slab[number];
        number++;
      }
    }
  } else{
    //cerr << "FluidcellTable<T>::set_xy_slab ERROR: size mismatch in received data, recv_slab.size()!=m_xcell*m_ycell - (" << recv_slab.size() << "!=" << m_xcell*m_ycell << ")"<< std::endl;
  }
}


/*!
  Get a value for each fluid cell using a fluid cell member function and return a vector
  of values, with position of cell.

  \param rdf the fluid cell member function
*/
template<typename T>
template<typename P>
vector<pair<Vec3,P> > NeighborTable<T>::forAllInnerCellsGet(P (CFluidCell::*rdf)() const)
{
  vector<pair<Vec3,P> > res;
  for(int ix=1;ix<m_xcell-1;ix++){
    for(int iy=1;iy<m_ycell-1;iy++){
      for(int iz=1;iz<m_zcell-1;iz++){
        int idx=cellindex(ix,iy,iz);
        Vec3 global_pos=m_cells[idx].getPos();
        res.push_back(make_pair(global_pos,(m_cells[idx].*rdf)()));
      }
    }
  }
  return res;
}


/*!
  Get a value for each fluid cell using a fluid cell member function and return a vector
  of values,with global index of cell

  \param rdf the fluid cell member function
*/
template<typename T>
template<typename P>
vector<pair<Vec3,P> > NeighborTable<T>::forAllInnerCellsGetIndexed(P (CFluidCell::*rdf)() const)
{
  vector<pair<Vec3,P> > res;
  for(int ix=1;ix<m_xcell-1;ix++){
    for(int iy=1;iy<m_ycell-1;iy++){
      for(int iz=1;iz<m_zcell-1;iz++){
        int idx=cellindex(ix,iy,iz);
        Vec3 global_index=m_cells[idx].getIndex();
        res.push_back(make_pair(global_index,(m_cells[idx].*rdf)()));
      }
    }
  }
  return res;
}


template<typename T>
template<typename P>
vector<P> NeighborTable<T>::forAllInnerCellsGetSum(P (CFluidCell::*rdf)() const)
{
  vector<P> res;
  for(int ix=1;ix<m_xcell-1;ix++){
    for(int iy=1;iy<m_ycell-1;iy++){
      for(int iz=1;iz<m_zcell-1;iz++){
        int idx=cellindex(ix,iy,iz);
        res.push_back((m_cells[idx].*rdf)());
      }
    }
  }
  return res;
}

/****fluid contents: end****/

