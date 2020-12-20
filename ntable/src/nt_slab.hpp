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
  Construct an empty NTSlab
*/ 
template<typename T>
NTSlab<T>::NTSlab():m_sl(0,0,0,0,0)
{
  m_table=NULL;
}

/*!
  Construct a NTSlab from a pointer to a NeigborTable
  and an valarray slice

  \param a the pointer to the search array
  \param s the  DSlice containing the slab parameters
*/
template<typename T>
NTSlab<T>::NTSlab(NeighborTable<T>* t,DSlice s):m_table(t),m_sl(s)
{
}
  
/*!
  get the number of particles in the slab
*/
template<typename T>
unsigned int NTSlab<T>::size()const
{
  unsigned int np=0;

  for(unsigned int i=0;i<m_sl.size();i++){
    np+=m_table->nparts_at_gridpoint(m_sl[i]);
  }

  return np;
}

/*!
  Return an iterator pointing to the first particle in the slab
  \todo check for end
*/
template<typename T>
typename NTSlab<T>::iterator NTSlab<T>::begin()
{
  // find the first occupied site in the search array
  unsigned int i=0;
  while(
    (i<m_sl.size())
    &&
    (m_table->nparts_at_gridpoint(m_sl[i]) == 0)
  ){
    i++;
  }
  return 
    typename NTSlab<T>::iterator(
      this,
      typename NeighborTable<T>::indextype(i,0)
    );
}

/*!
  Return an iterator pointing past the last particle in the slab
*/
template<typename T>
typename NTSlab<T>::iterator NTSlab<T>::end()
{
  return 
    typename NTSlab<T>::iterator(
      this,
      typename NeighborTable<T>::indextype(m_sl.size(),0)
    );
}

/*!
  Return an iterator pointing to the last particle in the slab
  \todo check for end
*/
template<typename T>
typename NTSlab<T>::iterator NTSlab<T>::rbegin()
{
  int j;

  // find the last occupied site in the search array
  int i=m_sl.size()-1;
  while((i>-1)&&
	nparts_at_gridpoint(i)==0){
    i--;
  }
  j=nparts_at_gridpoint(i)-1;
  if(i==-1) j=0; // rend

  return
    typename NTSlab<T>::iterator(
      this,
      typename NeighborTable<T>::indextype(i,j)
    ); 
}

/*!
  Return an iterator pointing before the first particle in the slab
*/
template<typename T>
typename NTSlab<T>::iterator NTSlab<T>::rend()
{
  return
    typename NTSlab<T>::iterator(
      this,
      typename NeighborTable<T>::indextype(-1,0)
    );
}

/*!
  Return pointer to particle at index. Does the translation from the 
  slab index to the array index.

  \param idx the index
*/
template<typename T>
T* NTSlab<T>::ptr(typename NeighborTable<T>::indextype idx)
{
  typename NeighborTable<T>::indextype a_idx;
  a_idx.first=m_sl[idx.first];
  a_idx.second=idx.second;

  return m_table->ptr(a_idx); 
}

/*!
  Return reference to particle at index. Does the translation from the 
  slab index to the array index.

  \param idx the index
*/
template<typename T>
T& NTSlab<T>::ref(typename NeighborTable<T>::indextype idx)
{ 
  typename NeighborTable<T>::indextype a_idx;
  a_idx.first=m_sl[idx.first];
  a_idx.second=idx.second;

  return m_table->ref(a_idx); 
}


/*!
  Insert single item. The iterator parameter is only there to ensure inteface compatibility with 
  STL containers and ignored otherwise.  

  \param data the item to be inserted
*/
template<typename T>
void NTSlab<T>::insert(typename NTSlab<T>::iterator,const T& data)
{
  m_table->insert(data);
}


/*!
  erase single item

  \param
*/
template<typename T>
void NTSlab<T>::erase(typename NTSlab<T>::iterator pos)
{
  // check if iter points into this slab ?
  // erase item 
  typename NeighborTable<T>::indextype a_idx;
  a_idx.first=m_sl[pos.index().first];
  a_idx.second=pos.index().second;
  m_table->erase(a_idx);
}

/*!
  erase range

  \param
  \param
*/
template<typename T>
void NTSlab<T>::erase(typename NTSlab<T>::iterator begin,typename NTSlab<T>::iterator end)
{
  // do it backwards, because erase invalidates all iterators behind
  begin--;
  end--;
  
  for(typename NTSlab<T>::iterator iter=end; 
      iter!=begin;
      iter--){
    
    erase(iter);
  }
}

template <typename T>
bool operator== (const NTSlab<T>& s1,const NTSlab<T>& s2)
{
  return(s1.m_table==s2.m_table && s1.m_sl==s2.m_sl);
}

template <typename T>
bool operator!= (const NTSlab<T>& s1,const NTSlab<T>& s2)
{
  return (s1.m_table!=s2.m_table || s1.m_sl!=s2.m_sl);
}

/*!
  output operator for NTSlab

  \param ost the output stream
  \param NTS the NTSlab
*/
template<typename T>
ostream& operator<< (ostream& ost,const NTSlab<T>& NTS)
{
  ost << "---NTSlab---" << endl;
  ost << "indices:" << endl;
  for(unsigned int i=0;i<NTS.slab_size();i++){
    ost << NTS.m_sl[i] << " ";
  }
  ost << endl;

  return ost;
}
