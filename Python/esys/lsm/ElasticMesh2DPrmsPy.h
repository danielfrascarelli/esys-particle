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
#ifndef ESYS_LSMELASTICMESH2DPRMSPY_H
#define ESYS_LSMELASTICMESH2DPRMSPY_H

#include "Model/ETriMeshIP.h"

#include <string>

using std::string; 

namespace esys
{
  namespace lsm
  {
    /*!
      \class NRotElasticMesh2DPrmsPy
      \brief class for elastic 2D mesh interactions in python interface. Deprecated: use NRotElasticLinMeshPrmsPy.
    */
    class NRotElasticMesh2DPrmsPy : public ETriMeshIP
    {
    public:
      NRotElasticMesh2DPrmsPy(
        const string &interactionName,
        const string &meshName,
        double        normalK
      );
    };
    
    /*!
      \class NRotElasticLinMeshPrmsPy
      \brief Class for elastic piece-wise linear mesh interactions in the python interface.
    */
    class NRotElasticLinMeshPrmsPy : public ETriMeshIP
    {
    public:
      NRotElasticLinMeshPrmsPy(
        const string &interactionName,
        const string &meshName,
        double        normalK
      );
    };
    
    void exportElasticMesh2DPrms();
  } // namespace lsm
} // namespace esys

#endif //ESYS_LSMELASTICMESH2DPRMSPY_H
