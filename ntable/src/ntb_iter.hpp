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
  Construct block iterator

  \param ntb the block over which to iterate
  \param ix the index in the x-dimension
  \param iy the index in the y-dimension
  \param iz the index in the z-dimension
  \param ig the index at the gridpoint
*/
template<typename T>
NTBlock_iter<T>::NTBlock_iter(NTBlock<T>* ntb,int ix,int iy,int iz,int ig)
  : m_block(ntb),
    m_ix(ix),
    m_iy(iy),
    m_iz(iz),
    m_ig(ig)
{
}

/*!
  prefix increment operator

  \warning not impl.
*/
template<typename T>
NTBlock_iter<T>& NTBlock_iter<T>::operator++()
{
}

/*!
  postfix increment operator
*/
template<typename T>
NTBlock_iter<T> NTBlock_iter<T>::operator++(int)
{
  NTBlock_iter<T> t=*this;

  if(m_ig+1<static_cast<int>(m_block->nparts_at_gridpoint(m_ix,m_iy,m_iz))){
    m_ig++;
  } else {
    m_ig=0;
    bool found=false;
    bool end=false;
    while(!found && !end){ 
      if(m_ix<m_block->m_xmax){
	m_ix++;
      } else {
	m_ix=m_block->m_xmin;
	if(m_iy<m_block->m_ymax){
	  m_iy++;
	} else {
	  m_iy=m_block->m_ymin;
	  if(m_iz<m_block->m_zmax){
	    m_iz++;
	  } else {
	    m_ix=m_block->m_xmax+1;// end iter 
	    m_iy=m_block->m_ymax+1;
	    m_iz=m_block->m_zmax+1;
	    end=true;
	  }
	}
      }
      if(!end) found=m_block->nparts_at_gridpoint(m_ix,m_iy,m_iz)>0;
    }
  }
  return t;
}

/*!
  prefix decrement operator

  \warning not impl.
*/
// template<typename T>
// NTBlock_iter<T>& NTBlock_iter<T>::operator--()
// {}

/*!
  postfix decrement operator
*/
// template<typename T>
// NTBlock_iter<T> NTBlock_iter<T>::operator--(int)
// {}

/*!
  access operator

  \todo what happens if end() is dereferenced ? 
*/
template<typename T>
T* NTBlock_iter<T>::operator->()
{
  return m_block->ptr(m_ix,m_iy,m_iz,m_ig);
}

/*!
  dereference  operator

  \todo what happens if end() is dereferenced ? 
*/
template<typename T>
T& NTBlock_iter<T>::operator*()
{
  return m_block->ref(m_ix,m_iy,m_iz,m_ig);
}

/*!
  equality operator
*/
template<typename T>
bool operator== (const NTBlock_iter<T>& b1,const NTBlock_iter<T>& b2)
{
  return (b1.m_block==b2.m_block &&
	  b1.m_ix==b2.m_ix &&
	  b1.m_iy==b2.m_iy &&
	  b1.m_iz==b2.m_iz &&
	  b1.m_ig==b2.m_ig);
}

/*!
  inequality operator
*/
template<typename T>
bool operator!= (const NTBlock_iter<T>& b1,const NTBlock_iter<T>& b2)
{ 
  return (b1.m_block!=b2.m_block ||
	  b1.m_ix!=b2.m_ix ||
	  b1.m_iy!=b2.m_iy ||
	  b1.m_iz!=b2.m_iz ||
	  b1.m_ig!=b2.m_ig);
}
