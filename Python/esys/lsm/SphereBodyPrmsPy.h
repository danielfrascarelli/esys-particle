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

#ifndef ESYS_LSMSPHEREBODYPRMSPY_H
#define ESYS_LSMSPHEREBODYPRMSPY_H

// --- project includes ---
#include "Python/esys/lsm/util/Vec3Py.h"
#include "Model/ESphereBodyInteractionGroup.h"



// --- STL includes ---
#include <string>

using namespace esys::lsm;

namespace esys
{
  namespace lsm
  {
    /*!
      \class NRotElasticSphereBodyPrmsPy
      \brief wrapper for CESphereBodyIGP 
      
      $Revision$
      $Date$
    */
    class NRotElasticSphereBodyPrmsPy : public  CESphereBodyIGP
    {
    public:
      NRotElasticSphereBodyPrmsPy(
        const std::string&,
        const std::string&,
        double
      );
    };

    void exportSphereBodyPrms();

  } // namespace lsm
} // namespace esys

#endif // ESYS_LSMSPHEREBODYPRMSPY_H
