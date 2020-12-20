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

#include <boost/version.hpp>
#include "Python/esys/lsm/TriggerPrmsPy.h"

namespace esys
{
  namespace lsm
  {
    MaxTriggerPrmsPy::MaxTriggerPrmsPy(double on,double off,int buff,int tail)
    {
      trig_on_value=on;
      trig_off_value=off;
      buff_size=buff;
      tail_size=tail;
    }

    using boost::python::arg;
    void exportTriggerPrms()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<MaxTriggerPrmsPy>(
	 "MaxTriggerPrms",
	 "Parameter class for max trigger",
	 boost::python::init<double,double,int,int>((
           arg("on_value"),
           arg("off_value"),
           arg("buffer_size"),
           arg("tail_size"))
         ));
	 					    
    }
  } // namespace lsm
} // namespace esys
