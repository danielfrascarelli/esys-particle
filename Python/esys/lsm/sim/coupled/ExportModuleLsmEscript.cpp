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


#include <boost/python.hpp>
#include "Python/esys/lsm/sim/coupled/LsmMpiEscriptPy.h"

BOOST_PYTHON_MODULE(LsmEscriptPy)
{
  esys::lsm::exportLsmMpiEscript();
};
