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


#ifndef ESYS_LSMTRIMESHBUILDPARAMSPY_H
#define ESYS_LSMTRIMESHBUILDPARAMSPY_H

/*---------------------------------------------------
 *
 * wrapper classes for TriMesh build parameters
 *
 *--------------------------------------------------*/

//--- project includes ---
#include "Model/BTriMeshIP.h"

namespace esys
{
  namespace lsm
  {
    /*!
      \class MeshTagBuildPrmsPy
      \brief wrapper for MeshTagBuildPrms
    */
    class MeshTagBuildPrmsPy : public MeshTagBuildPrms
    {
    public:
      MeshTagBuildPrmsPy(int, int);
    };

    /*!
      \class MeshGapBuildPrmsPy
      \brief wrapper for MeshGapBuildPrms
    */
    class MeshGapBuildPrmsPy : public MeshGapBuildPrms
    {
    public:
      MeshGapBuildPrmsPy(double);
    };

    void exportMeshBuildPrms();
  } // namespace lsm
} // namespace esys

#endif //ESYS_LSMTRIMESHBUILDPARAMSPY_H
