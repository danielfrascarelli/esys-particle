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

#ifndef __ERROR_H
#define __ERROR_H

#include <string>
#include <iostream>

using namespace std;

class MError
{
 protected:
  string message;

 public:
  MError(const string& m){message="Error : "+m;};
  virtual ~MError(){};

  void print_err() const{cerr <<  message << endl;};
};

#endif //__ERROR_H
