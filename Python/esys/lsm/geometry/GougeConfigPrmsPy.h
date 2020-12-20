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


#ifndef ESYS_LSM_GOUGEBLOCKPRMS_H
#define ESYS_LSM_GOUGEBLOCKPRMS_H

#include "Foundation/console.h"
#include "Python/esys/lsm/geometry/GougeConfigPy.h"
#include "Python/esys/lsm/util/BoundingBoxPy.h"

namespace boost
{
  namespace python
  {
    class list;
  }
}

namespace esys
{
  namespace lsm
  {
    typedef GougeConfigPy::GougeConfPrms GougeConfigPrmsBasePy;

    class ParticleRndPackPrmsPy : public GougeConfigPrmsBasePy::ParticleRndPackPrms
    {
    public:
      typedef GougeConfigPrmsBasePy::ParticleRndPackPrms Inherited;
      ParticleRndPackPrmsPy(double size, double minRadius, double maxRadius);
    };

    class GrainRndPackPrmsPy : public GougeConfigPrmsBasePy::GrainRPackPrms
    {
    public:
      typedef GougeConfigPrmsBasePy::GrainRPackPrms Inherited;
      typedef Inherited::ParticleGrainGen            ParticleGrainGenPy;
      
      GrainRndPackPrmsPy(
        double size,
        ParticleGrainGenPy &grainGen,
        int connTag
      );
    };

    class GougeConfigPrmsPy : public GougeConfigPrmsBasePy
    {
    public:
      typedef GougeConfigPrmsBasePy Inherited;
      GougeConfigPrmsPy(
        const BoundingBoxPy       &bBox,
        double                    padRadius,
        const ParticleRndPackPrmsPy        &roughnessPrms,
        const GrainRndPackPrmsPy   &gougePrms,
        int                       maxInsertionFailures,
        const boost::python::list &periodicDimList,
        double                    tolerance = DBL_EPSILON*128,
        double                    connectionTolerance = DBL_EPSILON*128*10,
        int                       blockConnTag=0
      );
    };

    void exportGougeConfigPrms();
  }
}

#endif
