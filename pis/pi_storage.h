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

#ifndef __PARALLEL_INTERACTION_STORAGE_H
#define __PARALLEL_INTERACTION_STORAGE_H

//--- MPI includes ---
#include <mpi.h>

//--- STL includes ---
#include <list>
#include <vector>
#include <iostream>
#include <stdexcept>

using std::list;
using std::vector;
using std::pair;

//--- Project includes ---
#include "Foundation/vec3.h"
#include "Foundation/quintuple.h"
#include "Foundation/triplet.h"
#include "Parallel/CheckPointable.h"
#include "tml/comm/comm.h"

// forward declared class to avoid circular include
class AFieldSlave; 
class AParallelParticleArray;

/*!
  \class AParallelInteractionStorage
  \brief abstract base class for parallel interaction storage array
*/
class AParallelInteractionStorage : public esys::lsm::CheckPointable
{
 protected:
  AParallelParticleArray* m_ppa;

 public:
  AParallelInteractionStorage(AParallelParticleArray* ppa){m_ppa=ppa;};
  virtual ~AParallelInteractionStorage(){};

  virtual void exchange()=0;
  virtual void rebuild()=0;
  virtual bool update()=0;
  //  virtual void tryInsert(const vector<int>&)=0;
  virtual bool isIn(const vector<int>&)=0;
  virtual void calcForces()=0;
  virtual void calcHeatFrict() {}
  virtual void calcHeatTrans() {}
  virtual void setTimeStepSize(double dt)=0;
  virtual void addExIG(AParallelInteractionStorage*){}; // do nothing
  virtual AFieldSlave* generateNewScalarFieldSlave(TML_Comm*,const string&,int,int,int,int)=0;
  virtual AFieldSlave* generateNewVectorFieldSlave(TML_Comm*,const string&,int,int,int,int)=0;
  virtual AFieldSlave* generateNewScalarHistoryFieldSlave(TML_Comm*,const string&,int,int,int);
  virtual void saveCheckPointData(std::ostream&)
  {
    throw std::runtime_error("saveCheckPointData not implemented in subclass.");
  }
  
  virtual void loadCheckPointData(std::istream&)
  {
    throw std::runtime_error("loadCheckPointData not implemented in subclass.");
  }

  virtual void saveSnapShotData(std::ostream&)
  {
    throw std::runtime_error(" saveSnapShotData not implemented in subclass.");
  }

  //!< access for setter functions
  virtual void forAllInteractionsSet(const string&,double){};

  virtual bool willSave(){ return false;};
};

/*!
  \class TParallelInteractionStorage
  \brief templated abstract base class for parallel interaction storage array.
  Adds the vector of interactions and access functions to AParallelInteractionStorage
*/
template <typename I>
class TParallelInteractionStorage : public AParallelInteractionStorage
{
 public:
  typedef I interaction_type;

 protected:
  list<I> m_interactions;

 public:
  TParallelInteractionStorage(AParallelParticleArray* ppa):AParallelInteractionStorage(ppa){};
  virtual ~TParallelInteractionStorage(){};

  class InteractionIterator {
  public:
    typedef I Interaction;
    typedef typename list<I>::iterator Iterator;

    InteractionIterator(Iterator begin, Iterator end, AParallelParticleArray* ppa);

    bool hasNext();

    Interaction &next();

    int getNumRemaining();
    
  protected:
    bool isInner(const Iterator &it);

  private:
    int                    m_numRemaining;
    Iterator               m_it;
    Iterator               m_end;
    AParallelParticleArray *m_ppa;
  };

  InteractionIterator getInnerInteractionIterator();

  //!< types
  typedef esys::lsm::quintuple<Vec3,double,Vec3,double,Vec3> Raw2Data;
  typedef esys::lsm::triplet<int,int,Vec3> DataWithID;
  typedef esys::lsm::quintuple<int,int,Vec3,Vec3,Vec3> DataWithPosID;

  //!< access functions
  template <typename P> vector<pair<Vec3,P> > forAllInnerInteractionsGetWithPos(P (I::*rdf)() const);
  template <typename P> vector<pair<Raw2Data,P> > forAllInnerInteractionsGetRaw2(P (I::*rdf)() const);
  template <typename P> vector<pair<DataWithID,P> > forAllInnerInteractionsGetDataWithID(P (I::*rdf)() const);
  template <typename P> vector<pair<DataWithPosID,P> > forAllInnerInteractionsGetDataWithPosID(P (I::*rdf)() const);
  template <typename P> void forAllInnerInteractionsGet(P&,typename P::value_type (I::*rdf)() const);

  //!< access for setter functions
  virtual void forAllInteractionsSet(const string&,double);

  //!< access functions with tags
  template <typename P> vector<pair<Vec3,P> > forAllTaggedInnerInteractionsGetWithPos(P (I::*rdf)() const,int,int);
  template <typename P> void forAllTaggedInnerInteractionsGet(P&,typename P::value_type (I::*rdf)() const,int,int);

  //!< generate FieldSlave of correct type
  virtual AFieldSlave* generateNewScalarFieldSlave(TML_Comm*,const string&,int,int,int,int); // const ?
  virtual AFieldSlave* generateNewVectorFieldSlave(TML_Comm*,const string&,int,int,int,int);
  
};

#include "pis/pi_storage.hpp"

#endif //__PARALLEL_INTERACTION_STORAGE_H
