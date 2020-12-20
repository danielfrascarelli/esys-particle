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
#include "Python/esys/lsm/ParticleIdPairVectorPy.h"
#include "Python/esys/lsm/util/BoundingBoxPy.h"
#include "Python/esys/lsm/geometry/SimpleSphereNeighboursPy.h"
#include "Python/esys/lsm/geometry/SimpleBlockPy.h"
#include "Python/esys/lsm/geometry/SimpleSpherePy.h"
#include "Python/esys/lsm/geometry/IteratorPy.h"
#include "Python/BoostPythonUtil/PythonIterIterator.h"
#include "Python/BoostPythonUtil/ListConverter.h"
#include "Geometry/SphereNeighbours.h"

namespace esys
{
  namespace lsm
  {
    class SimpleSphereNeighboursPy
      : public SphereNeighbours<SimpleSpherePy, ParticleIdPairVectorPy>
    {
    public:
      typedef
        SphereNeighbours<SimpleSpherePy, ParticleIdPairVectorPy>
        Inherited;

      SimpleSphereNeighboursPy(
        double maxDist,
        const BoundingBoxPy &bBox,
        const boost::python::list &circDimList
      )
        : Inherited(
            maxDist,
            bBox,
            bpu::listToVector<bool>(circDimList)
          )
      {
      }

      IdPairVector getNeighboursPy(boost::python::object &iteratable)
      {
        return
          getNeighbours(
            bpu::PythonIterIterator<SimpleSpherePy &>(iteratable)
          );
      }

      typedef IteratorPy<Inherited::Iterator> SSNIteratorPy;
      SSNIteratorPy getIteratorPy() const
      {
        return getIterator();
      }
    };

    using boost::python::arg;
    void exportSimpleSphereNeighbours()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<SimpleSphereNeighboursPy, boost::noncopyable>(
        "SimpleSphereNeighbours",
        "Discovers pairs of spheres which are closer than a specified\n"
        "threshold distance.\n",
        boost::python::init<
          double,
          const BoundingBoxPy &,
          const boost::python::list &
        >(
          (
            arg("maxDist") = 0.01,
            arg("bBox") = BoundingBoxPy(Vec3Py(-10,-10,-10), Vec3Py(10,10,10)),
            arg("circDimList") =
              boost::python::list(
                boost::python::make_tuple(false,false,false)
              )
          ),
          "Construct the neighbour-finding object.\n"
          "@type maxDist: float\n"
          "@kwarg maxDist: Threshhold distance governing which spheres"
          " are determined to be neighbours.\n"
          "@type bBox: L{BoundingBox<esys.lsm.util.FoundationPy.BoundingBox>}\n"
          "@kwarg bBox: Initial box for look-up grid. Also indicates locations"
          " of fixed circular boundaries.\n"
          "@type circDimList: list of 3 boolean elements\n"
          "@kwarg circDimList: Default tag given to the connections.\n"
        )
      )
      .def(
        "getNumSpheres",
        &SimpleSphereNeighboursPy::getNumSpheres,
        "Returns the number of spheres which have been added to this"
        " neighbour finder via the L{getNeighbours} method.\n"
        "@rtype: int\n"
        "@return: The number of spheres added to this neighbour finder."
      )
      .def(
        "getNumIdPairs",
        &SimpleSphereNeighboursPy::getNumIdPairs,
        "Returns the total number of neighbouring pairs of spheres "
        " found during calls to the L{getNeighbours} method.\n"
        "@rtype: int\n"
        "@return: number of C{ParticleIdPair} objects in this collection."
      )
      .def(
        "__iter__",
        &SimpleSphereNeighboursPy::getIteratorPy,
        "Iterator for enumerating sequence of"
        " all C{ParticleIdPair<esys.lsm.LsmPy.ParticleIdPair>}"
        " neighbours discovered during L{getNeighbours} calls.\n"
        "@rtype: L{SsNeighbourParticleIdPairIterator}\n"
        "@return: Iterator for enumerating sequence of"
        " C{ParticleIdPair<esys.lsm.LsmPy.ParticleIdPair>} objects."
      )
      .def("__len__", &SimpleSphereNeighboursPy::getNumIdPairs)
      .def(
        "getNeighbours",
        &SimpleSphereNeighboursPy::getNeighboursPy,
        (
          arg("iterable")
        ),
        "Returns a sequence of"
        " C{ParticleIdPair<esys.lsm.LsmPy.ParticleIdPair>} objects which"
        " indicate neighbouring pairs of particles. Only returns"
        " pairs of id's from particles within the specified sequence.\n"
        "@type iterable: iterable\n"
        "@kwarg iterable: Sequence of L{SimpleSphere} objects.\n"
        "@rtype: C{ParticleIdPairVector<esys.lsm.LsmPy.ParticleIdPairVector>}\n"
        "@return: Sequence of id pairs indicating spheres which are closer"
        " than the threshold distance.\n"
      )
      ;

      SimpleSphereNeighboursPy::SSNIteratorPy::exportIterator(
        "SsNeighbourParticleIdPairIterator"
      );
    }
  }
}
