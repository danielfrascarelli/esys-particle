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
#include <boost/python.hpp>
#include <boost/noncopyable.hpp>
#include "Foundation/console.h"
#include "Python/esys/lsm/geometry/PackerPy.h"
#include "Python/BoostPythonUtil/ListConverter.h"

namespace esys
{
  namespace lsm
  {
    //==========================================================================
    //==========================================================================
    //==========================================================================
    PackerPy::PackerPy(NTablePtr nTablePtr) : Inherited(nTablePtr)
    {
    }

    PackerPy::PackerPy(ParticlePoolPtr particlePoolPtr, NTablePtr nTablePtr)
      : Inherited(particlePoolPtr, nTablePtr)
    {
    }

    PackerPy::SimpleSphereIteratorPy PackerPy::getSimpleSphereIteratorPy()
    {
      return Inherited::getParticleIterator();
    }

    SimpleSphereCollectionPy PackerPy::getSimpleSphereCollectionPy()
    {
      return getParticleCollection();
    }

    void PackerPy::generate()
    {
      this->get_override("generate")();
    }

    //==========================================================================
    //==========================================================================
    //==========================================================================
    BoxPackerPy::BoxPackerPy(
      ParticlePoolPtr   particlePoolPtr,
      NTablePtr         nTablePtr,
      const BoundingBox &bBox,
      const BoolVector  &periodicDimensions,
      double            tolerance
    ) : Inherited(
          particlePoolPtr,
          nTablePtr,
          bBox,
          periodicDimensions,
          tolerance
        )
    {
    }
    //==========================================================================
    //==========================================================================
    //==========================================================================

    CubicBoxPackerPy::CubicBoxPackerPy(
      double              radius,
      const BoundingBoxPy &bBox,
      const boost::python::list &periodicDimensions,
      double              tolerance
    ) :
        Inherited(
          ParticleGeneratorPtr(
            new ParticleGenerator(radius)
          ),
          ParticlePoolPtr(new ParticlePool(4096)),
          NTablePtr(
            new NTable
            (
              bBox,
              (4.0*radius), // table grid size
              bpu::listToVector<bool>(periodicDimensions),
              4*radius
            )
          ),
          bBox,
          bpu::listToVector<bool>(periodicDimensions),
          tolerance,
          radius
        )
    {
    }
    //==========================================================================
    //==========================================================================
    //==========================================================================
    RandomSpherePackerPy::RandomSpherePackerPy(
      double minRadius,
      double maxRadius,
      double cubicPackRadius,
      int maxInsertionFailures,
      const BoundingSpherePy &bSphere,
      double tolerance,
      bool do2d
    )
      :
        Inherited(
          ParticleGeneratorPtr(new ParticleGenerator(minRadius, maxRadius)),
          ParticlePoolPtr(new ParticlePool(4096)),
          NTablePtr(
            new NTable
            (
              do2d ? bSphere.get2dBBox() : bSphere.getBBox(),
              (8.0*minRadius), // table grid size
              BoolVector(3, false),
              4*max(maxRadius, cubicPackRadius)
            )
          ),
          bSphere,
          tolerance,
          cubicPackRadius,
          maxInsertionFailures,
          do2d
        )
    {
    }

    //==========================================================================
    //==========================================================================
    //==========================================================================
    RandomBoxPackerPy::RandomBoxPackerPy(
      double minRadius,
      double maxRadius,
      double cubicPackRadius,
      int maxInsertionFailures,
      const BoundingBoxPy &bBox,
      const boost::python::list &periodicDimensions,
      double tolerance
    )
      :
        Inherited(
          ParticleGeneratorPtr(new ParticleGenerator(minRadius, maxRadius)),
          ParticlePoolPtr(new ParticlePool(4096)),
          NTablePtr(
            new NTable
            (
              bBox,
              (8.0*minRadius), // table grid size
              bpu::listToVector<bool>(periodicDimensions),
              4*max(maxRadius, cubicPackRadius)
            )
          ),
          bBox,
          bpu::listToVector<bool>(periodicDimensions),
          tolerance,
          cubicPackRadius,
          maxInsertionFailures
        )
    {
    }

    RandomBoxPackerPy::RandomBoxPackerPy(
      ParticleGeneratorPtr particleGeneratorPtr,
      ParticlePoolPtr      particlePoolPtr,
      NTablePtr            nTablePtr,
      const BoundingBox    &bBox,
      const BoolVector     &periodicDimensions,
      double               tolerance,
      double               cubicPackRadius,
      int                  maxInsertionFailures,
      const PlaneVector    &fitPlaneVector
    ) :
        Inherited(
          particleGeneratorPtr,
          particlePoolPtr,
          nTablePtr,
          bBox,
          periodicDimensions,
          tolerance,
          cubicPackRadius,
          maxInsertionFailures,
          fitPlaneVector
        )
    {
    }

    //==========================================================================
    //==========================================================================
    //==========================================================================
    RndGrainGenPy::RndGrainGenPy(
      double minGrainRadius,
      double maxGrainRadius,
      double minParticleRadius,
      double maxParticleRadius
    ) : Inherited(minGrainRadius, maxGrainRadius),
        m_minParticleRadius(minParticleRadius),
        m_maxParticleRadius(maxParticleRadius)
    {
    }

    const double &RndGrainGenPy::getMinParticleRadius() const
    {
      return m_minParticleRadius;
    }

    const double &RndGrainGenPy::getMaxParticleRadius() const
    {
      return m_maxParticleRadius;
    }

    RndGrainGenPy::Grain
    RndGrainGenPy::getGrain(const Particle &p)
    {
      return this->get_override("getGrain")(p);
    }

    GrainRandomBoxPackerPy::GrainRandomBoxPackerPy(
      ParticleGrainGen          &particleGrainGen,
      double                    cubicPackRadius,
      int                       maxInsertionFailures,
      const BoundingBox         &bBox,
      const boost::python::list &circDimList,
      double                    tolerance
    ) : Inherited(
          ParticleGrainGenPtr(),
          ParticlePoolPtr(new ParticlePool(4096)),
          NTablePtr(
            new NTable
            (
              bBox,
              (8.0*particleGrainGen.getMinGrainRadius()), // table grid size
              bpu::listToVector<bool>(circDimList),
              4*max(particleGrainGen.getMinGrainRadius(), cubicPackRadius)
            )
          ),
          bBox,
          bpu::listToVector<bool>(circDimList),
          tolerance,
          cubicPackRadius,
          maxInsertionFailures
        )
    {
      setParticleGrainGen(particleGrainGen);
    }

    GrainRandomBoxPackerPy::GrainRandomBoxPackerPy(
      ParticleGrainGenPtr particleGrainGenPtr,
      ParticlePoolPtr     particlePoolPtr,
      NTablePtr           nTablePtr,
      const BoundingBox   &bBox,
      const BoolVector    &periodicDimensions,
      double              tolerance,
      double              cubicPackRadius,
      int                 maxInsertionFailures,
      const PlaneVector   &fitPlaneVector,
      GrainPoolPtr        grainPoolPtr
    ) :
        Inherited(
          particleGrainGenPtr,
          particlePoolPtr,
          nTablePtr,
          bBox,
          periodicDimensions,
          tolerance,
          cubicPackRadius,
          maxInsertionFailures,
          fitPlaneVector,
          grainPoolPtr
        )
    {
    }

    
    GrainRandomBoxPackerPy::GrainIteratorPy
    GrainRandomBoxPackerPy::getGrainIteratorPy()
    {
      return Inherited::getGrainIterator();
    }

    const GrainRandomBoxPackerPy::GrainCollection &
    GrainRandomBoxPackerPy::getGrainCollectionPy() const
    {
      return Inherited::getGrainCollection();
    }
    
    //==========================================================================
    //==========================================================================
    //==========================================================================
    using boost::python::arg;
    using boost::python::args;
    using boost::python::return_value_policy;
    using boost::python::copy_const_reference;
    using boost::python::return_by_value;
    using boost::python::pure_virtual;
    void exportPacker()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<PackerPy,boost::noncopyable>(
        "Packer",
        "Base class for sphere packers."
        " The L{generate} method is over-ridden in subclasses to"
        " create sphere packings.",
        boost::python::no_init
      )
      .def(
        "getSimpleSphereIterator",
        &PackerPy::getSimpleSphereIteratorPy,
        "Returns an iterator which can be used to iterate over all"
        " spheres in the packing.\n"
        "@return: iterator"
        "@rtype: iterable\n"
      )
      .def(
        "getNumSpheres",
        &PackerPy::getNumParticles,
        "Returns the number of spheres generated in the packing.\n"
        "@rtype: int\n"
        "@return: number of spheres in packing."
      )
      .def(
        "getSimpleSphereCollection",
        &PackerPy::getSimpleSphereCollectionPy,
        "Returns generated spheres as a L{SimpleSphereCollection}.\n"
        "@rtype: L{SimpleSphereCollection}"
      )
      .def(
        "generate",
        boost::python::pure_virtual(&PackerPy::Inherited::generate),
        "Causes this I{packer} to create spheres in a packing."
      )
      ;
      boost::python::class_<
        BoxPackerBasePy,
        boost::python::bases<PackerPy>,
        boost::noncopyable
      >
      (
        "BoxPackerBase",
        "Base class for L{BoxPacker}.",
        boost::python::no_init
      );
      boost::python::class_<
        BoxPackerPy,
        boost::python::bases<BoxPackerBasePy>,
        boost::noncopyable
      >
      (
        "BoxPacker",
        "Base class of I{packers} which generate packings within a rectangular"
        " region.",
        boost::python::no_init
      )
      ;
      boost::python::class_<
        CubicBoxPackerBasePy,
        boost::python::bases<BoxPackerPy>,
        boost::noncopyable>
      (
        "CubicBoxPackerBase",
        "Base class for L{CubicBoxPacker}.",
        boost::python::no_init
      )
      ;

      boost::python::class_<
        CubicBoxPackerPy,
        boost::noncopyable,
        boost::python::bases<CubicBoxPackerBasePy>
      >(
        "CubicBoxPacker",
        "Instances generate a regular cubic close packing of identically"
        " sized spheres within the confines of a specified box.",
        boost::python::init<
          double,
          const BoundingBoxPy &,
          const boost::python::list &,
          double
        >(
          (
            arg("radius"),
            arg("bBox"),
            arg("circDimList")=boost::python::list(bpu::vectorToList(BoolVector(3,false))),
            arg("tolerance")=0.001
          )
          ,
          "Constructor.\n"
          "@type radius: float\n"
          "@kwarg radius: Radius of generated spheres.\n"
          "@type bBox: L{BoundingBox<esys.lsm.util.FoundationPy.BoundingBox>}\n"
          "@kwarg bBox: box specifying the region into which spheres"
          " are packed. A 2D packing can be generated by specifying a\n"
          " zero sized I{z} dimension, eg bBox=BoundingBox(Vec3(1,1,0),Vec3(21,21,0))\n"
          " will generate a 2D packing in the I{x}-I{y} plane.\n"
          "@type circDimList: list of 3 bool elements\n"
          "@kwarg circDimList: list indication which (if any) of the box dimensions\n"
          " is circular (note, only a single dimension may be circular).\n"
          "@type tolerance: float\n"
          "@kwarg tolerance: Generated spheres may overlap by no more than\n"
          " this amount.\n"
        )
      )
      .def(
        "generate",
        &CubicBoxPackerPy::generate,
        "Generates cubic packing of spheres."
      )
      ;
      //========================================================================
      boost::python::class_<
        RndCubicBoxPackerBasePy,
        boost::python::bases<BoxPackerPy>,
        boost::noncopyable>
      (
        "RndCubicBoxPackerBase",
        "Base class for L{RandomSpherePacker}.",
        boost::python::no_init
      )
      ;
      boost::python::class_<
        RandomSpherePackerBasePy,
        boost::python::bases<RndCubicBoxPackerBasePy>,
        boost::noncopyable>
      (
        "RandomSpherePackerBase",
        "Base class for L{RandomSpherePacker}.",
        boost::python::no_init
      )
      ;
      boost::python::class_<
        RandomSpherePackerPy,
        boost::noncopyable,
        boost::python::bases<RandomSpherePackerBasePy>
      >(
        "RandomSpherePacker",
        "Instances generate a packing of randomly sized spheres"
        " within the confines of a specified sphere.",
        boost::python::init<
          double,
          double,
          double,
          int,
          const BoundingSpherePy &,
          double,
          bool
        >(
          (
            arg("minRadius"),
            arg("maxRadius"),
            arg("cubicPackRadius"),
            arg("maxInsertFails"),
            arg("bSphere"),
            arg("tolerance")=0.001,
            arg("do2d")=false
          )
          ,
          "Constructor.\n"
          "@type minRadius: float\n"
          "@kwarg minRadius: Minimum radius of generated spheres.\n"
          "@type maxRadius: float\n"
          "@kwarg maxRadius: Maximum radius of generated spheres.\n"
          "@type cubicPackRadius: float\n"
          "@kwarg cubicPackRadius: Initial randomly sized I{seed} spheres\n"
          " are generated at positions corresponding to a regular cubic-packed grid\n"
          " of spheres with radius C{cubicPackRadius}.\n"
          "@type maxInsertFails: int\n"
          "@kwarg maxInsertFails: Stopping criterion for terminating the random.\n"
          " insertion algorithm. Algorithm terminates after C{maxInsertFails}\n"
          " number of attempts at fitting a randomly generated sphere into the\n"
          " box.\n"
          "@type bSphere: L{BoundingSphere<esys.lsm.util.FoundationPy.BoundingSphere>}\n"
          "@kwarg bSphere: Sphere specifying the region into which sub-spheres"
          " are packed.\n"
          "@type tolerance: float\n"
          "@kwarg tolerance: Generated spheres may overlap by no more than\n"
          " this amount.\n"
          "@type do2d: bool\n"
          "@kwarg do2d: All generated spheres are on the M{z=0} plane.\n"
        )
      )
      .def(
        "generate",
        &RandomSpherePackerPy::generate,
        "Generates packing of randomly sized spheres."
      )
      ;

      //========================================================================


      
      boost::python::class_<
        RandomBoxPackerBasePy,
        boost::python::bases<BoxPackerPy>,
        boost::noncopyable>
      (
        "RandomBoxPackerBase",
        "Base class for L{RandomBoxPacker}.",
        boost::python::no_init
      )
      ;
      boost::python::class_<
        RandomBoxPackerPy,
        boost::noncopyable,
        boost::python::bases<RandomBoxPackerBasePy>
      >(
        "RandomBoxPacker",
        "Instances generate a packing of randomly sized spheres"
        " within the confines of a specified box.",
        boost::python::init<
          double,
          double,
          double,
          int,
          const BoundingBoxPy &,
          const boost::python::list &,
          double
        >(
          (
            arg("minRadius"),
            arg("maxRadius"),
            arg("cubicPackRadius"),
            arg("maxInsertFails"),
            arg("bBox"),
            arg("circDimList")=boost::python::list(bpu::vectorToList(BoolVector(3,false))),
            arg("tolerance")=0.001
          )
          ,
          "Constructor.\n"
          "@type minRadius: float\n"
          "@kwarg minRadius: Minimum radius of generated spheres.\n"
          "@type maxRadius: float\n"
          "@kwarg maxRadius: Maximum radius of generated spheres.\n"
          "@type cubicPackRadius: float\n"
          "@kwarg cubicPackRadius: Initial randomly sized I{seed} spheres\n"
          " are generated at positions corresponding to a regular cubic-packed grid\n"
          " of spheres with radius C{cubicPackRadius}.\n"
          "@type maxInsertFails: int\n"
          "@kwarg maxInsertFails: Stopping criterion for terminating the random.\n"
          " insertion algorithm. Algorithm terminates after C{maxInsertFails}\n"
          " number of attempts at fitting a randomly generated spheres into the\n"
          " box.\n"
          "@type bBox: L{esys.lsm.util.FoundationPy.BoundingBox}\n"
          "@kwarg bBox: box specifying the region into which spheres"
          " are packed. A 2D packing can be generated by specifying a\n"
          " zero sized I{z} dimension, eg bBox=BoundingBox(Vec3(1,1,0),Vec3(21,21,0))\n"
          " will generate a 2D packing in the I{x}-I{y} plane.\n"
          "@type circDimList: list of 3 bool elements\n"
          "@kwarg circDimList: list indication which (if any) of the box dimensions\n"
          " is circular (note, only a single dimension may be circular).\n"
          "@type tolerance: float\n"
          "@kwarg tolerance: Generated spheres may overlap by no more than\n"
          " this amount.\n"
        )
      )
      .def(
        "generate",
        &RandomBoxPackerPy::generate,
        "Generates packing of randomly sized spheres."
      )
      ;

      boost::python::class_<RndGrainGenPy, boost::noncopyable>
      (
        "RndGrainGen",
        "Virtual class used to generate grains for L{GrainRandomBoxPacker} instances.",
        boost::python::init<double,double,double,double>
        (
          (
            arg("minGrainRadius"),
            arg("maxGrainRadius"),
            arg("minParticleRadius"),
            arg("maxParticleRadius")
          ),
          "Constructs a grain generator.\n"
          "@type minGrainRadius: float\n"
          "@kwarg minGrainRadius: minimum radius for a grain.\n"
          "@type maxGrainRadius: float\n"
          "@kwarg maxGrainRadius: maximum radius for a grain.\n"
          "@type minParticleRadius: float\n"
          "@kwarg minParticleRadius: minimum radius of a grain subsphere.\n"
          "@type maxParticleRadius: float\n"
          "@kwarg maxParticleRadius: maximum radius of a grain subsphere.\n"
        )
      )
      .def(
        "getMinGrainRadius",
        &RndGrainGenPy::getMinGrainRadius,
        return_value_policy<copy_const_reference>(),
        "Returns the desired minimum radius for a grain."
      )
      .def(
        "getMaxGrainRadius",
        &RndGrainGenPy::getMaxGrainRadius,
        return_value_policy<copy_const_reference>(),
        "Returns the desired maximum radius for a grain."
      )
      .def(
        "getMinParticleRadius",
        &RndGrainGenPy::getMinParticleRadius,
        return_value_policy<copy_const_reference>(),
        "Returns the minimum radius for grain subsphere."
      )
      .def(
        "getMaxParticleRadius",
        &RndGrainGenPy::getMaxParticleRadius,
        return_value_policy<copy_const_reference>(),
        "Returns the maximum radius for grain subsphere."
      )
      .def(
        "getGrain",
        pure_virtual(&RndGrainGenPy::Inherited::getGrain),
        (
          arg("sphere")
        )
      )
      // This second "getGrain" binding is only for storing the Epytext, so that
      // Epydoc will not report mark-up errors (which happens if the strings
      // are with the first "getGrain" (pure_virtual) binding).
      .def(
        "getGrain",
        &RndGrainGenPy::Inherited::getGrain,
        "Method to be overridden in child classes. This method is expected to\n"
        "return a L{Grain} object which lies entirely within a I{sphere}.\n"
        "The sphere object passed to this method has C{getRadius} and\n"
        "C{getCentre} methods.\n"
        "@type sphere: L{SimpleSphere}\n"
        "@kwarg sphere: Specifes position and maximum radius of grain which is\n"
        "to be returned by this method.\n"
        "@rtype: L{Grain}\n"
        "@return: A L{Grain} which lies entirely within the specified C{sphere}.\n"
      )
      ;
      
      boost::python::class_<
        GrainCubicBoxPackerPy,
        boost::python::bases<BoxPackerPy>,
        boost::noncopyable>
      (
        "GrainCubicBoxPacker",
        "Base class for L{GrainRndBoxPacker}.",
        boost::python::no_init
      )
      ;
      boost::python::class_<
        GrainRndBoxPackerPy,
        boost::python::bases<GrainCubicBoxPackerPy>,
        boost::noncopyable>
      (
        "GrainRndBoxPacker",
        "Base class for L{GrainRandomBoxPackerBase}.",
        boost::python::no_init
      )
      ;
      boost::python::class_<
        GrainRandomBoxPackerBasePy,
        boost::python::bases<GrainRndBoxPackerPy>,
        boost::noncopyable>
      (
        "GrainRandomBoxPackerBase",
        "Base class for L{GrainRandomBoxPacker}.",
        boost::python::no_init
      )
      ;
      boost::python::class_<
        GrainRandomBoxPackerPy,
        boost::noncopyable,
        boost::python::bases<GrainRandomBoxPackerBasePy>
      >(
        "GrainRandomBoxPacker",
        "Instances generate a packing of randomly sized aggregate grains"
        " within the confines of a specified box.",
        boost::python::init<
          RndGrainGenPy &,
          double,
          int,
          const BoundingBoxPy &,
          const boost::python::list &,
          double
        >(
          (
            arg("grainGenerator"),
            arg("cubicPackRadius"),
            arg("maxInsertFails"),
            arg("bBox"),
            arg("circDimList")=boost::python::list(bpu::vectorToList(BoolVector(3,false))),
            arg("tolerance")=0.001
          )
          ,
          "Constructor.\n"
          "@type grainGenerator: L{RndGrainGen}\n"
          "@kwarg grainGenerator: object for generating grains within a\n"
          " specifed sphere.\n"
          "@type cubicPackRadius: float\n"
          "@kwarg cubicPackRadius: Initial randomly sized I{seed} grains\n"
          " are generated at positions corresponding to a regular cubic-packed grid\n"
          " of spheres with radius C{cubicPackRadius}.\n"
          "@type maxInsertFails: int\n"
          "@kwarg maxInsertFails: Stopping criterion for terminating the random.\n"
          " insertion algorithm. Algorithm terminates after C{maxInsertFails}\n"
          " number of attempts at fitting a randomly generated sphere into the\n"
          " box.\n"
          "@type bBox: L{BoundingBox<esys.lsm.util.FoundationPy.BoundingBox>}\n"
          "@kwarg bBox: box specifying the region into which grains"
          " are packed. A 2D packing can be generated by specifying a\n"
          " zero sized I{z} dimension, eg bBox=BoundingBox(Vec3(1,1,0),Vec3(21,21,0))\n"
          " will generate a 2D packing in the I{x}-I{y} plane.\n"
          "@type circDimList: list of 3 bool elements\n"
          "@kwarg circDimList: list indicating which (if any) of the box dimensions\n"
          " is circular (note, only a single dimension may be circular).\n"
          "@type tolerance: float\n"
          "@kwarg tolerance: Generated spheres may overlap by no more than\n"
          " this amount.\n"
        )
      )
      .def(
        "generate",
        &GrainRandomBoxPackerPy::generate,
        "Generates packing of randomly sized grains."
      )
      .def(
        "getNumGrains",
        &GrainRandomBoxPackerPy::getNumGrains,
        "Returns the number of grains in the packing."
      )
      .def(
        "getGrainIterator",
        &GrainRandomBoxPackerPy::getGrainIteratorPy,
        "Returns an iterator for enumerating grains in the packing."
      )
      .def(
        "getGrainCollection",
        &GrainRandomBoxPackerPy::getGrainCollectionPy,
        "Returns the collection of grains in the packing.",
        return_value_policy<copy_const_reference>()
      )
      ;
    }
  }
}
