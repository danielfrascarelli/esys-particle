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
#include "Python/esys/lsm/geometry/SimpleBlockPy.h"
#include "Python/esys/lsm/geometry/SimpleSphereCollectionPy.h"
#include "Python/esys/lsm/geometry/SimpleSpherePy.h"
#include "Geometry/SimpleBlock.h"

namespace esys
{
  namespace lsm
  {
    SimpleBlockPy::SimpleBlockPy() : SimpleSphereCollectionPy()
    {
    }
    
    SimpleBlockPy::SimpleBlockPy(
      const boost::python::list &dimCount,
      double radius
    )
      : SimpleSphereCollectionPy()
    {
      SimpleBlockGenerator<SimpleSpherePy>
        generator = SimpleBlockGenerator<SimpleSpherePy>(
            boost::python::extract<int>(dimCount[0]),
            boost::python::extract<int>(dimCount[1]),
            boost::python::extract<int>(dimCount[2]),
            radius
        );
      generator.createParticles(*this);
    }

    using boost::python::arg;
    void exportSimpleBlock()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<
        SimpleBlockPy,
        boost::python::bases<SimpleSphereCollectionPy>,
        boost::noncopyable
      >(
        "SimpleBlock",
        "A regular packing of spheres.",
        boost::python::init<
          const boost::python::list &,
          boost::python::optional<double>
        >(
          (arg("dimCount"), arg("radius")),
          "Creates a collection of C{dimCount[0]*dimCount[2]*dimCount[2]}"
          " spheres arranged in a regular rectangular lattice packing.\n"
          "@type dimCount: list of 3 elements\n"
          "@kwarg dimCount: number of particles in each packing dimension\n"
          "@type radius: float\n"
          "@kwarg radius: Radius of spheres in the packing\n"
        )
      )
      .def(
        boost::python::init<>()
      )
      ;
    }
  }
}
