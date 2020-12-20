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

#ifndef __BASICINTERACTION_H
#define __BASICINTERACTION_H

//-- IO includes --
#include <iostream>

/*!
  \class BasicInteraction
  \brief Class to represent the common part of a pair 
  interaction, i.e. the IDs of the particles and the 
  interaction tag

  \author Steffen Abe
*/
class BasicInteraction
{
 public:
  typedef int Id;
  typedef int Tag;

 private:
  Id  m_p1;
  Id  m_p2;
  Tag m_tag;

 public:

  BasicInteraction(Id id1, Id id2,Tag tag=0);

  Id first() const {return m_p1;}
  
  Id second() const {return m_p2;}

  Id getP1Id() const
  {
    return first();
  }

  Id getP2Id() const
  {
    return second();
  }

  Tag getTag() const
  {
    return m_tag;
  }

  template <typename TmplVisitor>
  void visit(TmplVisitor &visitor) const
  {
    visitor.visitBasicInteraction(*this);
  }

  friend std::ostream& operator<<(std::ostream&,const BasicInteraction&);
  friend class BILess;
};

/*!
  \class BILess
  \brief function object for the ordering of BasicInteraction


  \author Steffen Abe
*/
class BILess
{
 public:
  bool operator()(const BasicInteraction&,const BasicInteraction&); 
};

#endif //__BASICINTERACTION_H
