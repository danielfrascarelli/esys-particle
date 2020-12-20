/////////////////////////////////////////////////////////////
//                                                         //
// Copyright (c) 2003-2014 by The University of Queensland //
// Centre for Geoscience Computing                         //
// http://earth.uq.edu.au/centre-geoscience-computing      //
//                                                         //
// Primary Business: Brisbane, Queensland, Australia       //
// Licensed under the Open Software License version 3.0    //
// http://www.apache.org/licenses/LICENSE-2.0          //
//                                                         //
/////////////////////////////////////////////////////////////


#ifndef ESYS_LSMFLUIDFIELDSAVERPRMSPY_H
#define ESYS_LSMFLUIDFIELDSAVERPRMSPY_H

#include "Python/esys/lsm/FieldSaverPrmsPy.h"

//--- STL includes ---
#include <string>

//--- Boost includes ---
#include <boost/python.hpp>

namespace esys
{
  namespace lsm
  {

    class FluidScalarFieldSaverPrmsPy : public FieldSaverPrmsPy
    {
    public:
      FluidScalarFieldSaverPrmsPy(
        const std::string &fieldName,
        const std::string &fileName,
        const std::string &fileFormat,
        int beginTimeStep,
        int endTimeStep,
        int timeStepIncr
      );
    };

    class FluidVectorFieldSaverPrmsPy : public FieldSaverPrmsPy
    {
    public:
      FluidVectorFieldSaverPrmsPy(
        const std::string &fieldName,
        const std::string &fileName,
        const std::string &fileFormat,
        int beginTimeStep,
        int endTimeStep,
        int timeStepIncr
      );
    };

    void exportFluidFieldSaverPrms();
  } // namespace lsm
} // namespace esys

#endif
