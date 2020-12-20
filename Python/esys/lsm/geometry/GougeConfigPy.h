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


#ifndef ESYS_LSM_GOUGEBLOCKPY_H
#define ESYS_LSM_GOUGEBLOCKPY_H

#include <boost/python.hpp>
#include "Foundation/console.h"
#include "Geometry/GougeConfig.h"
#include "Python/esys/lsm/geometry/PackerPy.h"
#include "Python/esys/lsm/util/BoundingBoxPy.h"
#include "Python/esys/lsm/geometry/TaggedIdConnectionPy.h"

namespace esys
{
  namespace lsm
  {
    class GougeConfigPrmsPy;
    class GougeConfigPy :
      public GougeConfig<GrainRandomBoxPackerPy,RandomBoxPackerPy,TaggedIdConnectionPy>
    {
    public:
      typedef
        GougeConfig<GrainRandomBoxPackerPy,RandomBoxPackerPy,TaggedIdConnectionPy>
        Inherited;

      GougeConfigPy(const GougeConfigPrmsPy &prms);

      void writeVtkXml(const std::string &fileName);

      boost::python::list getCircDimList() const;
      
      BoundingBoxPy getParticleBoundingBox();

      BoundingBoxPy getDomainBoundingBox();

      boost::python::list getConnectionList() const;

      class BBoxVisitor
      {
      public:
        BBoxVisitor();

        BoundingBoxPy getBBox() const;

        template <typename TmplParticle>
        void visitSimpleParticle(TmplParticle &particle);

      private:
        Vec3 m_minPt;
        Vec3 m_maxPt;
      };

    private:
    };

    void exportGougeConfig();
  }
}

#endif
