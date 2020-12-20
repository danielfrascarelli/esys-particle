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
  construct empty block
*/
template<typename T>
NTBlock<T>::NTBlock()
  : m_table(NULL),
    m_xmin(0),
    m_xmax(0),
    m_ymin(0),
    m_ymax(0),
    m_zmin(0),
    m_zmax(0)
{
}

/*!
  Construct a NTBlock from a pointer to a NeigborTable
  and the min and max inices in each dimension

  \param t the pointer to the NeighborTable
  \param xmin the minimum x index 
  \param xmax the maximum x index 
  \param ymin the minimum y index 
  \param ymax the maximum y index 
  \param zmin the minimum z index 
  \param zmax the maximum z index 
*/
template<typename T>
NTBlock<T>::NTBlock(NeighborTable<T>* t,int xmin,int xmax,int ymin,int ymax,int zmin,int zmax)
  : m_table(t),
    m_xmin(xmin),
    m_xmax(xmax),
    m_ymin(ymin),
    m_ymax(ymax),
    m_zmin(zmin),
    m_zmax(zmax)
{
}

/*!
  get the number of particles in the block
*/
template<typename T>
unsigned int NTBlock<T>::size()
{
  unsigned int np=0;
  for (iterator it = begin(); it != end(); it++, np++)
  {
    // keep going
  }

  return np;
}


/*!
  Return an iterator pointing to the first particle in the slab
  \todo check for end
*/
template<typename T>
typename NTBlock<T>::iterator NTBlock<T>::begin()
{
  // find the first occupied site in the search array
  int ix=m_xmin;
  int iy=m_ymin;
  int iz=m_zmin;
  bool end=false;
  while(!end && m_table->nparts_at_gridpoint(m_table->index(ix,iy,iz))==0){
    if(ix<m_xmax){
      ix++;
    } else {
      ix=m_xmin;
      if(iy<m_ymax){
	iy++;
      } else {
	iy=m_ymin;
	if(iz<m_zmax){
	  iz++;
	} else {
	  ix=m_xmax+1;// end iter 
	  iy=m_ymax+1;
	  iz=m_zmax+1;
	  end=true;
	}
      }
    }
  }
  return iterator(this,ix,iy,iz,0);
}

/*!
  Return an iterator pointing past the last particle in the slab
*/
template<typename T>
typename NTBlock<T>::iterator NTBlock<T>::end()
{
  return iterator(this,m_xmax+1,m_ymax+1,m_zmax+1,0);
}

/*!
  Return an iterator pointing to the last particle in the slab
  \todo check for end
*/
// template<typename T>
// iterator NTBlock<T>::rbegin()
// {
  
// }

/*!
  Return an iterator pointing before the first particle in the slab
*/
// template<typename T>
// iterator NTBlock<T>::rend()
// {}
  
/*!
  Return pointer to particle at index.
*/
template<typename T>
T* NTBlock<T>::ptr(int ix,int iy,int iz,int j)
{
  typename NeighborTable<T>::indextype a_idx;
  a_idx.first=m_table->index(ix,iy,iz);
  a_idx.second=j;

  return m_table->ptr(a_idx); 
}

/*!
  Return reference to particle at index.
*/
template<typename T>
T& NTBlock<T>::ref(int ix,int iy,int iz,int j)
{
  typename NeighborTable<T>::indextype a_idx;
  a_idx.first=m_table->index(ix,iy,iz);
  a_idx.second=j;

  return m_table->ref(a_idx); 
}

/*!
  equality operator
*/
template <typename T>
bool operator== (const NTBlock<T>& b1,const NTBlock<T>& b2)
{
  return (b1.m_table==b2.m_table &&
	  b1.m_xmin==b2.m_xmin &&
	  b1.m_xmax==b2.m_xmax &&
	  b1.m_ymin==b2.m_ymin &&
	  b1.m_ymax==b2.m_ymax &&
	  b1.m_zmin==b2.m_zmin &&
	  b1.m_zmax==b2.m_zmax);
}

/*!
  inequality operator
*/
template <typename T>
bool operator!=(const NTBlock<T>& b1,const NTBlock<T>& b2)
{
  return (b1.m_table!=b2.m_table ||
	  b1.m_xmin!=b2.m_xmin ||
	  b1.m_xmax!=b2.m_xmax ||
	  b1.m_ymin!=b2.m_ymin ||
	  b1.m_ymax!=b2.m_ymax ||
	  b1.m_zmin!=b2.m_zmin ||
	  b1.m_zmax!=b2.m_zmax);
}


/*!
  output operator for NTBlock

  \param ost the output stream
  \param NTB the NTBlock
*/
template <typename T>
ostream& operator<<(ostream& ost,const NTBlock<T>& NTB)
{
  ost << "---NTBlock---"<< endl;
  ost << "range: " << endl;
  ost << "x: " << NTB.m_xmin << "-" << NTB.m_xmax << endl;
  ost << "y: " << NTB.m_ymin << "-" << NTB.m_ymax << endl;
  ost << "z: " << NTB.m_zmin << "-" << NTB.m_zmax << endl;
  ost << "indices:" << endl;
  for (int ix=NTB.m_xmin;ix<=NTB.m_xmax;ix++){    
    for (int iy=NTB.m_ymin;iy<=NTB.m_ymax;iy++){
      for (int iz=NTB.m_zmin;iz<=NTB.m_zmax;iz++){
	ost << NTB.m_table->index(ix,iy,iz) <<  " ";
      }
    }
  }
  ost << endl;

  return ost;
}
