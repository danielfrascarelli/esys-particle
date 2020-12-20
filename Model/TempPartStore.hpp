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

/*!
  Construct a new TTempPartStore

  \param min minimum corner of the volume
  \param max maximum corner of the volume
  \param nx nr. of slaves in x-direction
  \param ny nr. of slaves in y-direction
  \param nz nr. of slaves in z-direction
*/
template<typename T>
TTempPartStore<T>::TTempPartStore(const Vec3& min,const Vec3& max, int nx,int ny,int nz)
{
  // minimum corner
  m_xmin=min.X();
  m_ymin=min.Y();
  m_zmin=min.Z();
  // size of a block
  m_xsize=(max.X()-min.X())/double(nx);
  m_ysize=(max.Y()-min.Y())/double(ny);
  m_zsize=(max.Z()-min.Z())/double(nz);
  // number of blocks
  m_nx=nx;
  m_ny=ny;
  m_nz=nz;
}

/*!
  helper function which translates block coordinates into an index

  \param x the x coordinate
  \param y the y coordinate
  \param z the z coordinate
  \return the index
*/
template<typename T>
int TTempPartStore<T>::coordsToIndex(int x,int y,int z)
{
  return x*m_ny*m_nz+y*m_nz+z;
}

/*!
  helper function which returns the index of the block in which a given point is located.

  \param pos the position of the point
  \return the index
*/
template<typename T>
int TTempPartStore<T>::posToIndex(const Vec3& pos)
{
  int posx=int(floor((pos.X()-m_xmin)/m_xsize));
  int posy=int(floor((pos.Y()-m_ymin)/m_ysize));
  int posz=int(floor((pos.Z()-m_zmin)/m_zsize));
  return coordsToIndex(posx,posy,posz);
}

/*!
  add a new slave to the coordinate->rank mapping table

  \param cx x-coordinate of the slave
  \param cy y-coordinate of the slave
  \param cz z-coordinate of the slave
  \param rank the rank of the slave (in MPI_COMM_WORLD)
*/
template<typename T>
void TTempPartStore<T>::addSlaveID(int cx,int cy,int cz,int rank)
{
  int idx=coordsToIndex(cx,cy,cz);
  m_slave_id_map.insert(make_pair(idx,rank));
}

/*!
  add a new particle to the storage

  \param p the particle
*/
template<typename T>
void TTempPartStore<T>::addParticle(const T& p)
{
  int idx=posToIndex(p.getPos());
  int sl_id=m_slave_id_map[idx];
  cout << "Inserting particle at pos " << p.getPos() << " index " << idx << "slave " << sl_id << endl;
  typename multimap<int,T>::iterator it=m_mmap.insert(make_pair(sl_id,p));
  m_by_id.insert(make_pair(p.getID(),it));
}

/*!
  add a connection between2 particles to the storage

  \param id1 the Id of the first particle
  \param id2 the Id of the second particle
  \param tag the connection tag

  \warning not implemented
*/
template<typename T>
void TTempPartStore<T>::addConnection(int id1,int id2,int tag)
{}
