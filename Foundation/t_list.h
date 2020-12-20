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

#ifndef _T_LIST_H_
#define _T_LIST_H_
#include <string.h>
#ifndef NULL
#define NULL 0
#endif

#ifdef _LINUX
#include <stddef.h>
#endif

#ifdef _DEBUG
#include "console.h"
#endif

template<class T>
class Node
{
public:
  Node<T> *Next, *Prev ;
  T *Val ;
} ;

/*!
  List container.
*/
template <class T>
class List
   {
   protected:
     Node<T> *Start ;   //!< Pointer to Start of list
     Node<T> *End;      //!< Pointer to end of list
     Node<T> *Current ; //! Pointer to current position in the list
   public:
      inline List() ;
      inline List(const List &L) ;
      virtual ~List() ;

      inline void Swap() ;
      inline void InsertAtStart(T *V) ;	
      inline void Append(T *V) ;	
      inline void InsertAfter(T *V) ; 
      inline void InsertBefore(T *V) ;  
      inline T* Get() ;
      inline void Put(T *V) ;

      inline void Clear() ;
      inline void Destroy() ;

      inline List& operator << (T *V) ;
      inline List& operator >> (T *V) ;

      inline void Next() ;
      inline void Prev() ;
      inline void First() ;
      inline void Last() ;

      inline int IsEnd() ;
      inline int IsStart() ;

      inline int SizeList() ;

      inline List operator + (const List& L) ;	
      inline List& operator += (const List& L) ;
      inline List& operator = (const List& L) ;
   } ;


/*!
  Stack container.
*/
template< class T>
class Stack {
    protected:
	List<T> L ;
    public:
	inline Stack() {} ;
	virtual ~Stack() { L.Destroy() ; } ;

	inline void Push(T* V) ;

	inline T* Pop() ;
} ;

#if 0
/*!
  List container where current position in list can be save
  using pop and push.
*/
template <class T>
class ListWS: public List<T> {
    protected:
	Stack<Node<T> > stack ;
    public:
	inline ListWS() : List<T>() {} ;
	inline ListWS(const ListWS &L) : List<T>(L) { } ;
	virtual ~ListWS() { } ;

	inline void PushPos() ; //!< Save current position in the list
	inline void PopPos() ; //!< restore last saved position in the list
} ;
#endif
#include "t_list.hpp"

#endif 

