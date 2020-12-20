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

#ifndef __PARALLEL_PARTICLE_ARRAY_H
#define __PARALLEL_PARTICLE_ARRAY_H

//--- MPI ---
#include <mpi.h>

//--- project includes ---
#include "ntable/src/ntable.h"
#include "ntable/src/nt_block.h"
#include "tml/comm/comm.h"
#include "tml/comm/cart_comm.h"
#include "Foundation/vec3.h"
#include "Geometry/Triangle.h"
#include "Geometry/AEdge.h"
#include "Foundation/quadtuple.h"
#include "Model/FluidCell.h" //fluid contents
#include "tml/comm/cart_comm.h" //fluid contents

//--- STL includes ---
#include <vector>
#include <set>

using std::vector;
using std::set;

//--- IO includes ---

/*!
  \class AParallelParticleArray
  \brief abstract base class for parallel particle storage array
 */
class AParallelParticleArray
{
 protected:
  TML_CartComm m_comm;
  int m_timestamp;

 public:
  AParallelParticleArray(TML_Comm *comm, const std::vector<unsigned int> &dims);
  AParallelParticleArray(TML_Comm *comm, const std::vector<unsigned int> &dims, const std::vector<bool> &circ);
  // virtual destructor
  virtual ~AParallelParticleArray(){};

  // get communicator
  TML_CartComm getComm() const {return m_comm;};

  //! return time stamp of last rebuild
  int getTimeStamp(){return m_timestamp;};

  // get ids of boundary particles
  virtual set<int> getBoundarySlabIds(int,int) const=0;
  virtual set<int> get2ndSlabIds(int,int) const=0;

  // check if pos is in inner part
  virtual bool isInInner(const Vec3&)=0;
};


/*!
  \class ParallelParticleArray
  \brief parrallel particle storage array with neighborsearch and
  variable exchange
*/
template<typename T>
class ParallelParticleArray : public AParallelParticleArray
{
 public: // types
  typedef T_Handle<typename NeighborTable<T>::pairlist> PairListHandle;
  typedef typename NeighborTable<T>::pairlist::iterator PairListIterator;
  typedef T_Handle<typename NeighborTable<T>::particlelist> ParticleListHandle;
  typedef typename NeighborTable<T>::particlelist::iterator ParticleListIterator;

  typedef T_Handle<typename NeighborTable<T>::particlecelllist> ParticleCellListHandle; //fluid contents
  typedef typename NeighborTable<T>::particlecelllist::iterator ParticleCellListIterator;//fluid contents

 private:
  NeighborTable<T>* m_nt;
  Vec3 m_minpos,m_maxpos; //!< local minimum and maximum positions
  double m_xshift,m_yshift,m_zshift; //!< circular shift values
  bool m_circ_edge_x_up,m_circ_edge_x_down; //!< circular edge flags
  static const int m_exchg_tag;
  bool m_fluidinitiated; //fluid contents
  bool m_firststep; //fluid contents

  // helper fnc
  template<typename P> void exchange_single(P (T::*rdf)(),void (T::*wrtf)(const P&),NTSlab<T>,NTSlab<T>,int,int);
  
 protected:

 public:
  ParallelParticleArray(TML_Comm *comm, const vector<unsigned int> &dims,const Vec3 &min,const Vec3 &max, double rmax,double alpha);
  ParallelParticleArray(TML_Comm *comm, const vector<unsigned int> &dims, const vector<bool> &circ,const Vec3 &min,const Vec3 &max, double rmax,double alpha);
  ~ParallelParticleArray();

  // info func
  Vec3 getMinPos()const {return m_minpos;};
  Vec3 getMaxPos()const {return m_maxpos;};
  vector<int> getCommCoords() const {return m_comm.get_coords();};
  vector<int> getCommDims() const {return m_comm.get_all_dims();};
  int size(){return m_nt->size();};
  int getInnerSize(){return (m_nt->inner()).size();};
  list<esys::lsm::quadtuple<int,int,int,int> > getOccupiedCellPositions(){return m_nt->getOccupiedCellPositions();};

  /****fluid contents: begin***/
  void setFluid(double,double,double,double,double,double,double,Vec3,Vec3);
  void setFluidVec3(Vec3,double,double,double,double,double,double,Vec3,Vec3);
  bool fluidInitiated(){return m_fluidinitiated;};
  void updateFluid();
  void setPressure(vector<pair<Vec3,double> >);
  void exchangeCells(double,int);
  void exchange_boundary(vector<CFluidCell>, int, int);
  template<typename P> vector<pair<Vec3,P> > forAllInnerCellsGet(P (CFluidCell::*rdf)() const);
  template<typename P> vector<pair<Vec3,P> > forAllInnerCellsGetIndexed(P (CFluidCell::*rdf)() const);
  template<typename P> vector<P> forAllInnerCellsGetSum(P (CFluidCell::*rdf)() const);
  /****fluid contents: end***/

  // particle insert ops
  void insert(const T&); //!< particle insertion
  void insert(const vector<T>&); //!< multi particle insert

  // check if pos is in inner part
  virtual bool isInInner(const Vec3&);

  // particle access (ugly!)
  T* getParticlePtrByIndex(int);
  T* getParticlePtrByPosition(const Vec3&);

  // rebuild
  void rebuild();

  //--- collective particle ops ---
  // variable exchange
  template<typename P> void exchange(P (T::*rdf)(),void (T::*wrtf)(const P&));

  // call member func for single particle by id
  void forParticle(int,void (T::*rdf)());
  template <typename P> void forParticle(int,void (T::*rdf)(P),const P&);

  // call member func for single particle by tag
  void forParticleTag(int,void (T::*rdf)());
  template <typename P> void forParticleTag(int,void (T::*rdf)(P),const P&);
  void forParticleTagMask(int,int,void (T::*rdf)());
  template <typename P> void forParticleTagMask(int,int,void (T::*rdf)(P),const P&);

  // call member func for all particles, different nr. of params
  void forAllParticles(void (T::*rdf)());
  void forAllParticles(void (T::*rdf)()const);
  template <typename P> void forAllParticles(void (T::*rdf)(P),const P&);

  // call member func for all inner particles
  template <typename P> void forAllInnerParticles(void (T::*rdf)(P&),P&);

  class ParticleIterator
  {
  public:
    typedef NTBlock<T> NtBlock;
    typedef T Particle;
    typedef typename NtBlock::iterator BlockIterator;

    ParticleIterator(const NtBlock &ntBlock);

    bool hasNext() const;

    Particle &next();

    int getNumRemaining() const;

  private:
    NtBlock       m_ntBlock;
    BlockIterator m_it;
    int           m_numRemaining;
  };

  ParticleIterator getInnerParticleIterator();

  // particle data access functions
  template <typename P> void forAllParticlesGet(P&,typename P::value_type (T::*rdf)() const);
  template <typename P> void forAllInnerParticlesGet(P&,typename P::value_type (T::*rdf)() const);
  template <typename P> vector<pair<int,P> > forAllParticlesGetIndexed(P (T::*rdf)() const);
  template <typename P> vector<pair<int,P> > forAllInnerParticlesGetIndexed(P (T::*rdf)() const);

  // particle data access functions with tag check
  template <typename P> void forAllTaggedParticlesGet(P&,typename P::value_type (T::*rdf)() const,int,int);
  template <typename P> void forAllTaggedInnerParticlesGet(P&,typename P::value_type (T::*rdf)() const,int,int);
  template <typename P> vector<pair<int,P> > forAllTaggedParticlesGetIndexed(P (T::*rdf)() const,int,int);
  template <typename P> vector<pair<int,P> > forAllInnerTaggedParticlesGetIndexed(P (T::*rdf)() const,int,int);

  // geometric data access function
  template <typename P> void forPointsGetNearest(P&,typename P::value_type (T::*rdf)() const,const Vec3&,double,double,double,int,int,int);

  // get ids of boundary particles
  virtual set<int> getBoundarySlabIds(int,int) const;
  virtual set<int> get2ndSlabIds(int,int) const;

  //--- get neigborlist stuff ---
  //! Get list of all pairs. Forwards to NTable::getFullList().
  PairListHandle getFullPairList(){return m_nt->getFullList();};
  //! Get list of new pairs. Forwards to NTable::getNewList().
  PairListHandle getNewPairList(){return m_nt->getNewList();};
  //! Get list of particles along a plane. Forwards to NTable::getParticlesAtPlane
  ParticleListHandle getParticlesAtPlane(Vec3 o,Vec3 n){return m_nt->getParticlesAtPlane(o,n);};
  //! Get list of particles near a sphere body. Forwards to NTable::getParticlesNearSphere
  ParticleListHandle getParticlesNearSphere(Vec3 c,double r){return m_nt->getParticlesNearSphere(c,r);};
  //! Get list of particles near a triangle. Forwards to NTable::getParticlesNearTriangle
  ParticleListHandle getParticlesNearTriangle(const Triangle& t){return m_nt->getParticlesNearTriangle(t);};
  //! Get list of particles near an edge. Forwards to NTable::getParticlesNearEdge
  ParticleListHandle getParticlesNearEdge(const AEdge* e){return m_nt->getParticlesNearEdge(e);};
  //! Get list of particles near a point. Forwards to NTable::getParticlesNearEdge
  ParticleListHandle getParticlesNearPoint(const Vec3& v){return m_nt->getParticlesNearPoint(v);};
  //! Get list of all particles. Forwards to NTable
  ParticleListHandle getAllParticles(){return m_nt->getAllParticles();};

  //! get all particles in inner block and put them into a vector
  void getAllInnerParticles(vector<T>&);

  /****fluid contents: begin***/
  //! Get list of all pairs of particle and cell. Forwards to NTable::getParticleCellList().
  ParticleCellListHandle getParticleCellPairList(){return m_nt->getParticleCellList();};
  //! Get list of new pairs of particle and cell. Forwards to NTable::getNewParticleCellList().
  ParticleCellListHandle getNewParticleCellPairList(){return m_nt->getNewParticleCellList();};
  /****fluid contents: end***/

  //--- checkpointing ---
  void saveCheckPointData(std::ostream&);
  void loadCheckPointData(std::istream&);

  //--- output (for debugging)---
  template <typename TT>
  friend ostream& operator<<(ostream &, const ParallelParticleArray<TT> &);
};

#include "ppa/src/pp_array.hpp"

#endif //__PARALLEL_PARTICLE_ARRAY_H

