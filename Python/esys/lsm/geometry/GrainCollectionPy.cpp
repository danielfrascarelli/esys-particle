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
#include "Python/esys/lsm/geometry/GrainCollectionPy.h"
#include "Python/esys/lsm/geometry/SimpleSpherePy.h"

namespace esys
{
  namespace lsm
  {
    GrainCollectionPy::GrainCollectionPy()
      : Inherited()
    {
    }

    GrainCollectionPy::GrainCollectionPy(ParticlePoolPtr particlePoolPtr)
      : Inherited(particlePoolPtr)
    {
    }

    GrainCollectionPy::GrainCollectionPy(
      ParticlePoolPtr particlePoolPtr,
      GrainPoolPtr    grainPoolPtr
    )
      : Inherited(particlePoolPtr, grainPoolPtr)
    {
    }


    GrainPy &GrainCollectionPy::createGrainPy()
    {
      return Inherited::createGrain();
    }

    GrainPy &GrainCollectionPy::createGrainWithIdPy(GrainPy::Id id)
    {
      return Inherited::createGrain(id);
    }

    GrainCollectionPy::GrainIteratorPy GrainCollectionPy::getGrainIteratorPy()
    {
      return GrainIteratorPy(Inherited::getGrainIterator());
    }

    void exportGrainCollection()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<GrainCollectionPy>(
        "GrainCollection",
        boost::python::init<>()
      )
      .def("getGrainIterator", &GrainCollectionPy::getGrainIteratorPy)
      .def("getNumParticles", &GrainCollectionPy::getNumParticles)
      .def("getNumGrains", &GrainCollectionPy::getNumGrains)
      .def("__iter__", &GrainCollectionPy::getGrainIteratorPy)
      .def("__len__", &GrainCollectionPy::getNumGrains)
      .def(
        "createGrain",
        &GrainCollectionPy::createGrainPy,
        boost::python::return_value_policy<
          boost::python::reference_existing_object
        >()
      )
      .def(
        "createGrain",
        &GrainCollectionPy::createGrainWithIdPy,
        boost::python::return_value_policy<
          boost::python::reference_existing_object
        >()
      )
      ;

      GrainCollectionPy::GrainIteratorPy::exportIterator(
        "GrainCollectionGrainIterator"
      );
    }
  }
}
