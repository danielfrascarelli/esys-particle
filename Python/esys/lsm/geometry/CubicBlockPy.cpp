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
#include "Python/esys/lsm/geometry/CubicBlockPy.h"
#include "Python/esys/lsm/geometry/SimpleSphereCollectionPy.h"
#include "Python/esys/lsm/geometry/SimpleSpherePy.h"
#include "Geometry/CubicBlock.h"

namespace esys
{
  namespace lsm
  {
    CubicBlockPy::CubicBlockPy() : SimpleSphereCollectionPy()
    {
    }
    
    CubicBlockPy::CubicBlockPy(
      const boost::python::list &dimCount,
      double radius,
      const ClosePackOrientation &orientation
    )
      : SimpleSphereCollectionPy()
    {
      CubicBlock<SimpleSpherePy>::BlockGenerator
        generator = CubicBlock<SimpleSpherePy>::BlockGenerator(
            boost::python::extract<int>(dimCount[0]),
            boost::python::extract<int>(dimCount[1]),
            boost::python::extract<int>(dimCount[2]),
            radius,
            orientation
        );
      generator.createParticles(*this);
    }

    using boost::python::arg;
    void exportOrientation()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::enum_<ClosePackOrientation>(
        "Orientation"
      )
        .value("DEFAULT", DEFAULT_ORIENT)
        .value("XYZ", XYZ)
        .value("XZY", XZY)
        .value("YXZ", YXZ)
        .value("YZX", YZX)
        .value("ZXY", ZXY)
        .value("ZYX", ZYX)
      ;
    }
    
    void exportCubicBlock()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<
        CubicBlockPy,
        boost::python::bases<SimpleSphereCollectionPy>,
        boost::noncopyable
      >(
        "CubicBlock",
        "Collection of particles arranged in a cubic close "
        " (face centred cubic) packing.",
        boost::python::init<
          const boost::python::list &,
          double,
          const ClosePackOrientation &
        >(
          (
            arg("dimCount"),
            arg("radius")=0.5,
            arg("orientation")=DEFAULT_ORIENT
          )
        )
      )
      .def(
        boost::python::init<
        >(
          "Creates a collection of C{dimCount[0]*dimCount[2]*dimCount[2]}"
          " spheres arranged in a cubic close packing.\n"
          "@type dimCount: list of 3 elements\n"
          "@kwarg dimCount: number of particles in each packing dimension\n"
          "@type radius: float\n"
          "@kwarg radius: Radius of spheres in the packing\n"
          "@type orientation: L{Orientation}\n"
          "@kwarg orientation: orientation of the axis-aligned layers"
          " in the packing."
        )
      )
      ;
    }
  }
}
