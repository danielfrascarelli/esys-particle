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

#ifndef __MESH2D_PIS_EB_H
#define __MESH2D_PIS_EB_H

#include "tml/comm/cart_comm.h"
#include "pis/mesh2d_pis.h"
#include <vector>

/*!
  \brief Class for parallel storage of interactions between a 2D mesh
  and particles which does require exchange of interactions across
  process boundaries but where interactions are not dynamically formed
*/
template<class ParticleType,class IType>
class Mesh2D_PIS_EB  : public Mesh2D_PIS<ParticleType>
{
 private:
  static const int m_edge_exchg_tag;
  static const int m_corner_exchg_tag;
  void exchange_boundary(int,int);

 protected:
  typename IType::ParameterType m_param;

  TML_CartComm m_comm;
  std::set<pair<int,int> > m_edge_int_set; // for isIn, <edge ID,particle ID> pairs 
  std::set<pair<int,int> > m_corner_int_set; // for isIn, <corner ID,particle ID> pairs 

  std::list<typename IType::TriIntType> m_edge_interactions;
  std::list<typename IType::CornerIntType> m_corner_interactions;

 public:

  // --- iterator class ===
  class InteractionIterator {
  public:
    typedef typename IType::TriIntType Interaction;
    typedef typename list<Interaction>::iterator Iterator;

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

  // --- member functions ---
  Mesh2D_PIS_EB(Mesh2D*,ParallelParticleArray<ParticleType>*,typename IType::ParameterType);
   
  virtual bool isIn(const vector<int>&);
  virtual void setTimeStepSize(double dt);
  virtual void calcForces();
  virtual bool update();
  virtual void exchange();
  virtual void rebuild();
  virtual void tryInsert(const typename IType::TriIntType&);
  virtual void tryInsert(const typename IType::CornerIntType&);
  virtual void tryInsert(const vector<int>&);

  InteractionIterator getInnerInteractionIterator();

  void buildFromPPATagged(int,int);
  void buildFromPPAByGap(double);
  
  virtual void saveSnapShotData(std::ostream&);
  virtual void saveCheckPointData(std::ostream&);
  virtual void loadCheckPointData(std::istream&);
};

#include "pis/mesh2d_pis_eb.hpp"

#endif //__MESH2D_PIS_EB_H
