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
#include "Python/esys/lsm/geometry/GrainPy.h"
#include "Python/esys/lsm/geometry/SimpleSpherePy.h"

namespace esys
{
  namespace lsm
  {
    GrainPy::GrainPy()
      : Inherited()
    {
    }

    GrainPy::GrainPy(Id id)
      : Inherited(id)
    {
    }

    GrainPy::GrainPy(ParticlePoolPtr particlePoolPtr)
      : Inherited(particlePoolPtr)
    {
    }

    GrainPy::GrainPy(Id id, ParticlePoolPtr particlePoolPtr)
      : Inherited(id, particlePoolPtr)
    {
    }

    boost::python::tuple
    GrainPy::getinitargs()
    {
      return boost::python::make_tuple(getId());
    }

    using boost::python::arg;
    void exportGrain()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<
        GrainPy,
        boost::python::bases<ParticleCollectionPy>
      >(
        "Grain",
        "Represents an aggregation of particles.",
        boost::python::init<>(
          "Constructs an empty grain with default id."
        )
      )
      .def(
        boost::python::init<GrainPy::Id>(
          ( boost::python::arg("id") ),
          "Constructs an empty grain with specified Id.\n"
          "@type id: int\n"
          "@kwarg id: Grain Id."
        )
      )
      .def(
        "getId",
        &GrainPy::getId,
        "Returns the Id of this grain.\n"
        "@rtype: int\n"
        "@return: Grain Id."
      )
      .def(
        "setId",
        &GrainPy::setId,
        ( arg("id") ),
        "Sets the id of this grain.\n"
        "@type id: int\n"
        "@kwarg id: The new grain Id."
      )
      .def("__getinitargs__", &GrainPy::getinitargs)
      ;
    }
  }
}
