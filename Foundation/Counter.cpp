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

#include "Counter.h"
#include "console.h"
#include <string.h>

CCounter::CCounter() {
    Name = NULL ;
    Value = 0 ;
    Id = 0 ;
}
CCounter::CCounter(char *iName, int iId) {
    create(iName,iId) ;
}
CCounter::~CCounter() {
    if (Name) delete Name ;
}

void CCounter::create(char *iName, int iId) {
    Name = new char[strlen(iName)+1] ;
    strcpy(Name,iName) ;
    Id = iId ;
    Value = 0 ;
}

CCounter & CCounter::operator +=(int n) {
    Value += n ;
    return *this ;
}
CCounter & CCounter::operator -=(int n) {
    Value -= n ;
    return *this ;
}
CCounter & CCounter::operator ++() {
    Value += 1 ;
    return *this ;
}
CCounter & CCounter::operator --() {
    Value -= 1 ;
    return *this ;
}
CCounter & CCounter::operator ++(int) {
    Value += 1 ;
    return *this ;
}
CCounter & CCounter::operator --(int) {
    Value -= 1 ;
    return *this ;
}
CCounter & CCounter::reset() {
    Value = 0 ;
    return *this ;
}

char *CCounter::getName() {
    return Name ;
}
    
CCounter::operator int () {
    return Value ;
}

CCounterList::CCounterList() {
}

CCounterList::~CCounterList() {
    m_Counters.First() ;
    while (!m_Counters.IsEnd()) {
        delete m_Counters.Get() ;
        m_Counters.Clear() ;
    }
}

CCounterList & CCounterList::operator << (CCounter &Counter) {
    CCounter *newC = new CCounter(Counter.getName()) ;
    m_Counters.Append(newC) ;
    return *this ;
}

CCounter & CCounterList::counter(char *name) {
    m_Counters.First() ;
    while (!m_Counters.IsEnd()) {
        if (!strcmp(name,m_Counters.Get()->getName())) {
            return *m_Counters.Get() ;
        }
        m_Counters.Next() ;
    }
    console.Error() << "Internal Error, Cannot find Counter " << name << "\n" ;
    CCounter *newC = new CCounter ;
    return *newC ;
}

CCounterList & CCounterList::addCounter(char *name) {
    CCounter *newC = new CCounter(name) ;
    m_Counters.Append(newC) ;
    return *this ;
}

CCounter & CCounterList::operator()(char *name) {
    return counter(name) ;
}

ostream&  CCounterList::print(ostream& Out) {
    m_Counters.First() ;
    while (!m_Counters.IsEnd()) {
        m_Counters.Get()->print(Out) ;
        Out << "\n" ;
        m_Counters.Next() ;
    }
    return Out ;
}

ostream&  CCounter::print(ostream& Out) {
    Out << Name << " = " << Value ;
    return Out ;
}

ostream&  operator<<(ostream& Out, CCounter &P) {
    return P.print(Out) ;
}

ostream&  operator<<(ostream& Out, CCounterList &P) {
    return P.print(Out) ;
}
