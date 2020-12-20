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

#ifndef _COUNTER_H_
#define _COUNTER_H_
#include "t_list.h"

//--- IO includes ---
#include <iostream>
#include <cstdio>
using std::ostream;

//--- system includes ---
#include <string.h>

/*!
  Provide a basic "counter", a counter has a Name and
  a value and an optional Id.
*/
class CCounter {
protected:
    char *Name ; //!< Name is allocated or deallocated on construction or destruction, respectively
    int Value ; 
    int Id ;
public:
    CCounter() ;
    CCounter(char *Name, int Id=0) ;
    virtual ~CCounter() ;

    void create(char *Name, int Id=0) ; //!< create counter, call by constructor
    CCounter & operator +=(int n) ;     //!< increment counter by n
    CCounter & operator -=(int n) ;     //!< decrement counter by n
    CCounter & operator ++() ;          //!< increment counter by 1
    CCounter & operator --() ;          //!< decrement counter by 1
    CCounter & operator ++(int) ;       //!< increment counter by 1
    CCounter & operator --(int) ;       //!< decrement counter by 1
    CCounter & reset() ;                //!< reset the counter to 0
    char *getName() ;                   //!< return name of counter
    operator int () ;                   //!< return the value of the counter
    ostream&  print(ostream& Out) ;     //!< method to print the value 
} ;
/*!
  List of CCounter (avoid duplication of code).
 
*/
class CListCounters : public List<CCounter> {
} ;

/*!
  Provide a list of counters.
*/
class CCounterList {
protected:
    CListCounters m_Counters ;
public:
    CCounterList() ;
    virtual ~CCounterList() ;

    CCounterList & operator << (CCounter &Counter) ;//!< add a counter "Counter"
    CCounterList & addCounter(char *name) ;         //!< add a counter of name "name"
    CCounter & counter(char *name) ;                //!< return the counter "name"
    CCounter & operator()(char *Name) ;             //!< return the counter "name"
    ostream&  print(ostream& Out) ;                 //!< method to print out all counters
    inline CListCounters & getList()                //!< return the list of counters
        { return m_Counters; } ; 
} ;

// out-of-class method to print-out counter values
ostream&  operator<<(ostream& Out, CCounter &P) ;
ostream&  operator<<(ostream& Out, CCounterList &P) ;

#endif

