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

#ifndef __IGPARAM_H
#define __IGPARAM_H

#include "Parallel/mpivbuf.h"

// -- STL includes --
#include <string>

/*!
  \class AIGParam
  \brief Abstract base class for InteractionGroup parameters
  \author Steffen Abe
  $Revision$
  $Date$
 */
class AIGParam
{
private:
  std::string m_name;

public:
  AIGParam(const std::string &name = "");
  
  virtual ~AIGParam();
  
  virtual void packInto(CVarMPIBuffer*) const;

  void setName(const std::string &name);
  
  const std::string &getName() const {return m_name;}
  
  const std::string &Name() const {return getName();}

  virtual std::string getTypeString() const = 0;

  //  friend class AInteractionGroup;
};

#endif //__IGPARAM_H
