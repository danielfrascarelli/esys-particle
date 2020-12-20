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

#ifndef __INTERACTIONPARAM_H
#define __INTERACTIONPARAM_H

// -- STL includes --
#include <string>
#include <map>

using std::string;
using std::map;

/*!
  \class AIParam
  \brief Abstract base class for interaction parameters 

  \author Steffen Abe
  $Revision$
  $Date$
*/
class AIParam
{
private:
  string m_name;
  map<string,double> m_data;

public:
  AIParam(const string&);
  virtual ~AIParam()=0;

  string getName()const{return m_name;};
  double getParamByName(const string&);
  void addParameter(const string&,double);
};

#endif //__INTERACTIONPARAM_H
