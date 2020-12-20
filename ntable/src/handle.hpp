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

template<typename T>
T_Handle<T>::T_Handle(T* t): m_rep(t), m_count(new int(1))
{}

template<typename T>
T_Handle<T>::T_Handle(const T_Handle& h): m_rep(h.m_rep), m_count(h.m_count)
{
   (*m_count)++;
}

template<typename T>
T_Handle<T>::~T_Handle()
{
  if(--(*m_count)==0){
    delete m_rep;
    delete m_count;
  }
}

template<typename T>
T_Handle<T>& T_Handle<T>::operator=(const T_Handle& h)
{
  if(m_rep!=h.m_rep){
    if(--(*m_count)==0){
      delete m_rep;
      delete m_count;
    }
    m_rep=h.m_rep;
    m_count=h.m_count;
    (*m_count)++;
  }
  return *this;
}

/*template<typename T>
void T_Handle<T>::destroy()
{
  if (m_count!=1) {
    throw HandleException;
  } else {
    
  }
  }*/
