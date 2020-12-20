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

#ifndef ESYS_LSMPACKERPY_H
#define ESYS_LSMPACKERPY_H

#include "Geometry/Packer.h"
#include "Geometry/BoxPacker.h"
#include "Geometry/CubicBoxPacker.h"
#include "Geometry/RandomBoxPacker.h"
#include "Geometry/RandomSpherePacker.h"
#include "Geometry/GrainRandomBoxPacker.h"
#include "Geometry/PackerGenerators.h"
#include "Python/esys/lsm/util/BoundingBoxPy.h"
#include "Python/esys/lsm/util/BoundingSpherePy.h"
#include "Python/esys/lsm/util/Vec3Py.h"
#include "Python/esys/lsm/geometry/IteratorPy.h"
#include "Python/esys/lsm/geometry/SimpleSphereCollectionPy.h"
#include "Python/esys/lsm/geometry/GrainPy.h"
#include "Python/esys/lsm/geometry/GrainCollectionPy.h"

namespace boost
{
  namespace python
  {
    class list;
    class tuple;
  }
}

namespace esys
{
  namespace lsm
  {
    class PackerPy
      : public Packer<SimpleSphereCollectionPy>,
        public boost::python::wrapper<Packer<SimpleSphereCollectionPy> >
    {
    public:
      typedef Packer<SimpleSphereCollectionPy> Inherited;
      typedef
        SimpleSphereCollectionPy::SimpleSphereIteratorPy
        SimpleSphereIteratorPy;

      PackerPy(NTablePtr nTablePtr);

      PackerPy(ParticlePoolPtr particlePoolPtr, NTablePtr nTablePtr);

      SimpleSphereIteratorPy getSimpleSphereIteratorPy();

      SimpleSphereCollectionPy getSimpleSphereCollectionPy();

      virtual void generate();
    };

    typedef BoxPacker<PackerPy> BoxPackerBasePy;

    class BoxPackerPy : public BoxPackerBasePy
    {
    public:
      typedef BoxPackerBasePy Inherited;

      BoxPackerPy(
        ParticlePoolPtr   particlePoolPtr,
        NTablePtr         nTablePtr,
        const BoundingBox &bBox,
        const BoolVector  &periodicDimensions,
        double            tolerance
      );
    };

    typedef ConstRadiusGen<SimpleSpherePy>              ConstRadiusGenPy;
    typedef CubicBoxPacker<ConstRadiusGenPy, BoxPackerPy> CubicBoxPackerBasePy;
    class CubicBoxPackerPy : public CubicBoxPackerBasePy
    {
    public:
      typedef CubicBoxPackerBasePy Inherited;
      CubicBoxPackerPy(
        double radius,
        const BoundingBoxPy &bBox,
        const boost::python::list &periodicDimensions,
        double              tolerance
      );
    };

    typedef RndRadiusGen<SimpleSpherePy>   RndRadiusGenPy;
    template <typename TPartGen>
    class PackerWrap
    {
    public:
      typedef CubicBoxPacker<TPartGen, BoxPackerPy> CubicBoxPackerBase;
      typedef RandomBoxPacker<TPartGen, ::esys::lsm::PackerWrap> RandomBoxPackerBase;
    };
    
    typedef PackerWrap<RndRadiusGenPy>::CubicBoxPackerBase  RndCubicBoxPackerBasePy;
    typedef PackerWrap<RndRadiusGenPy>::RandomBoxPackerBase RandomBoxPackerBasePy;
    class RandomBoxPackerPy : public RandomBoxPackerBasePy
    {
    public:
      typedef RandomBoxPackerBasePy Inherited;
      
      RandomBoxPackerPy(
        double minRadius,
        double maxRadius,
        double cubicPackRadius,
        int maxInsertionFailures,
        const BoundingBoxPy &bBox,
        const boost::python::list &periodicDimensions,
        double tolerance
      );

      RandomBoxPackerPy(
        ParticleGeneratorPtr particleGeneratorPtr,
        ParticlePoolPtr      particlePoolPtr,
        NTablePtr            nTablePtr,
        const BoundingBox    &bBox,
        const BoolVector     &periodicDimensions,
        double               tolerance,
        double               cubicPackRadius,
        int                  maxInsertionFailures,
        const PlaneVector    &fitPlaneVector
      );
    };

    typedef RandomSpherePacker<RndRadiusGenPy, PackerWrap> RandomSpherePackerBasePy;
    class RandomSpherePackerPy : public RandomSpherePackerBasePy
    {
    public:
      typedef RandomSpherePackerBasePy Inherited;
      
      RandomSpherePackerPy(
        double minRadius,
        double maxRadius,
        double cubicPackRadius,
        int maxInsertionFailures,
        const BoundingSpherePy &bSphere,
        double tolerance,
        bool do2d
      );
    };

    typedef GrainRndRadiusGen<GrainPy> RndGrainGenBasePy;
    class RndGrainGenPy :
      public RndGrainGenBasePy,
      public boost::python::wrapper<RndGrainGenBasePy>
    {
    public:
      typedef RndGrainGenBasePy Inherited;

      RndGrainGenPy(
        double minGrainRadius,
        double maxGrainRadius,
        double minParticleRadius,
        double maxParticleRadius
      );

      const double &getMinParticleRadius() const;

      const double &getMaxParticleRadius() const;

      virtual Grain getGrain(const Particle &p);

    private:
      double m_minParticleRadius;
      double m_maxParticleRadius;
    };

    typedef PackerWrap<RndGrainGenPy>::CubicBoxPackerBase  GrainCubicBoxPackerPy;
    typedef PackerWrap<RndGrainGenPy>::RandomBoxPackerBase GrainRndBoxPackerPy;
    typedef
      GrainRandomBoxPacker<RndGrainGenPy,GrainCollectionPy,PackerWrap>
      GrainRandomBoxPackerBasePy;
    class GrainRandomBoxPackerPy : public GrainRandomBoxPackerBasePy
    {
    public:
      typedef GrainRandomBoxPackerBasePy Inherited;
      typedef Inherited::Grain           Grain;
      typedef
        Inherited::GrainCollection::GrainIteratorPy
        GrainIteratorPy;

      GrainRandomBoxPackerPy(
        ParticleGrainGen          &particleGrainGen,
        double                    cubicPackRadius,
        int                       maxInsertionFailures,
        const BoundingBox         &bBox,
        const boost::python::list &circDimList,
        double                    tolerance
      );

      GrainRandomBoxPackerPy(
        ParticleGrainGenPtr particleGrainGenPtr,
        ParticlePoolPtr     particlePoolPtr,
        NTablePtr           nTablePtr,
        const BoundingBox   &bBox,
        const BoolVector    &periodicDimensions,
        double              tolerance,
        double              cubicPackRadius,
        int                 maxInsertionFailures,
        const PlaneVector   &fitPlaneVector,
        GrainPoolPtr         grainPoolPtr
      );


      GrainIteratorPy getGrainIteratorPy();

      const GrainCollection &getGrainCollectionPy() const;
    };

    void exportPacker();
  }
}

#endif
