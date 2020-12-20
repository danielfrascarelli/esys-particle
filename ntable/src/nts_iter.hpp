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
  Construct a NTSlab iterator
*/
template <typename T>
NTSlab_iter<T>::NTSlab_iter(NTSlab<T> *a,typename NeighborTable<T>::indextype idx):
  m_slab(a),m_curr(idx)
{}

/*!
  prefix increment operator

  \warning not impl.
*/
template <typename T>
NTSlab_iter<T>& NTSlab_iter<T>::operator++()
{
  return *this;
}

/*!
  postfix increment operator
*/
template <typename T>
NTSlab_iter<T> NTSlab_iter<T>::operator++(int)
{
  NTSlab_iter<T> t=*this;

  if(m_curr.second+1 < static_cast<int>(m_slab->nparts_at_gridpoint(m_curr.first))) {
    m_curr.second++;
  } else {
    m_curr.second=0;
    m_curr.first++;
    while(
      (m_curr.first < static_cast<int>(m_slab->slab_size()))
      &&
	    (m_slab->nparts_at_gridpoint(m_curr.first) == 0)
    ){
      m_curr.first++;
    }
  }

  return t;
}

/*!
  prefix decrement operator

  \warning not impl.
*/
template <typename T>
NTSlab_iter<T>& NTSlab_iter<T>::operator--()
{
  return *this;
}

/*!
  postfix decrement operator
*/
template <typename T>
NTSlab_iter<T> NTSlab_iter<T>::operator--(int)
{
  NTSlab_iter<T> t=*this;

  if(m_curr.second>0) {
    m_curr.second--;
  } else {
    m_curr.first--;
    while((m_curr.first>=0)&&
	  m_slab->nparts_at_gridpoint(m_curr.first)==0){
      m_curr.first--;
    }
    if(m_curr.first!=-1){
      m_curr.second=m_slab->nparts_at_gridpoint(m_curr.first)-1;
    }  
  }

  return t;
}

/*!
  access operator

  \todo what happens if end() is dereferenced ? 
*/
template <typename T>
T* NTSlab_iter<T>::operator->()
{
  return m_slab->ptr(m_curr);
}

/*!
  dereference  operator

  \todo what happens if end() is dereferenced ? 
*/
template <typename T>
T& NTSlab_iter<T>::operator*()
{
  return m_slab->ref(m_curr);
}

/* 
   get current index
*/
template <typename T>
typename NeighborTable<T>::indextype NTSlab_iter<T>::index() const
{
  return m_curr;
}

/*!
  equality operator
*/
template <typename T>
bool operator== (const NTSlab_iter<T>& i1,const NTSlab_iter<T>& i2)
{
  return (i1.m_slab==i2.m_slab && 
	  i1.m_curr.first==i2.m_curr.first &&
	  i1.m_curr.second==i2.m_curr.second);
}

/*!
  inequality operator
*/
template <typename T>
bool operator!= (const NTSlab_iter<T>& i1,const NTSlab_iter<T>& i2)
{
  return (i1.m_slab!=i2.m_slab || 
	  i1.m_curr.first!=i2.m_curr.first ||
	  i1.m_curr.second!=i2.m_curr.second);

}
