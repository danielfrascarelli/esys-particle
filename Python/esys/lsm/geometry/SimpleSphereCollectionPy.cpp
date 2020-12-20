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
#include "Python/esys/lsm/geometry/SimpleSphereCollectionPy.h"
#include "Python/esys/lsm/geometry/SimpleSpherePy.h"

namespace esys
{
  namespace lsm
  {
    SimpleSphereCollectionPy::SimpleSphereCollectionPy()
      : Inherited()
    {
    }

    SimpleSphereCollectionPy::SimpleSphereCollectionPy(
      const Inherited &particleCollection
    )
      : Inherited(particleCollection)
    {
    }

    SimpleSphereCollectionPy::SimpleSphereCollectionPy(
      const SimpleSphereCollectionPy &ssphereCollection
    )
      : Inherited(ssphereCollection)
    {
    }

    SimpleSphereCollectionPy::SimpleSphereCollectionPy(
      ParticlePoolPtr particlePoolPtr
    )
      : Inherited(particlePoolPtr)
    {
    }

    BoundingBoxPy SimpleSphereCollectionPy::getParticleBBoxPy() const
    {
      return
        BoundingBoxPy(
          ParticleCollection<SimpleSpherePy>::getParticleBBox()
        );
    }

    SimpleSphereCollectionPy::SimpleSphereIteratorPy
    SimpleSphereCollectionPy::getSimpleSphereIteratorPy()
    {
      return
        SimpleSphereIteratorPy(
          ParticleCollection<SimpleSpherePy>::getParticleIterator()
        );
    }

    SimpleSpherePy &SimpleSphereCollectionPy::createParticlePy(const SimpleSpherePy &p)
    {
      return ParticleCollection<SimpleSpherePy>::createParticle(p);
    }

    void SimpleSphereCollectionPy::rotatePy(const Vec3Py &rotation, const Vec3Py &pt)
    {
      rotate(rotation, pt);
    }

    void SimpleSphereCollectionPy::translateByPy(const Vec3Py &translation)
    {
      translateBy(translation);
    }

    class SimpleSphereCollectionPyPickleSuite : public boost::python::pickle_suite
    {
    public:
      static
      boost::python::tuple
      getstate(boost::python::object pcObj)
      {
        const SimpleSphereCollectionPy &pc =
          boost::python::extract<const SimpleSphereCollectionPy &>(pcObj);
        boost::python::list l;
        SimpleSphereCollectionPy::ParticleConstIterator it = pc.getParticleIterator();
        while (it.hasNext())
        {
          l.append(it.next());
        }
        return boost::python::make_tuple(pcObj.attr("__dict__"), l);
      }

      static
      void
      setstate(boost::python::object pcObj, boost::python::tuple state)
      {
        // restore the object's __dict__
        boost::python::dict d =
          boost::python::extract<boost::python::dict>(pcObj.attr("__dict__"))();
        d.update(state[0]);

        SimpleSphereCollectionPy &pc =
          boost::python::extract<SimpleSphereCollectionPy &>(pcObj);
        boost::python::list l = boost::python::extract<boost::python::list>(state[1]);
        for (int i = 0; i < bpu::len(l); i++)
        {
          pc.createParticle(
            boost::python::extract<SimpleSphereCollectionPy::Particle &>(l[i])
          );
        }
      }

      static bool getstate_manages_dict()
      {
        return true;
      }
    };

    using boost::python::arg;
    void exportSimpleSphereCollection()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<SimpleSphereCollectionPy>(
        "SimpleSphereCollection",
        "A collection of L{SimpleSphere} objects.",
        boost::python::init<>()
      )
      .def(
        "getBBox",
        &SimpleSphereCollectionPy::getParticleBBoxPy,
        "Returns the axis aligned bounding box which contains all spheres in\n"
        "this collection.\n"
        "@rtype: L{esys.lsm.util.BoundingBox<esys.lsm.util.FoundationPy.BoundingBox>}\n"
        "@return: tight bounding box for all spheres in this collection."
      )
      .def(
        "getNumSpheres",
        &SimpleSphereCollectionPy::getNumParticles,
        "Returns the number of spheres in this collection.\n"
        "@rtype: int\n"
        "@return: number of spheres in this collection"
      )
      .def(
        "__iter__",
        &SimpleSphereCollectionPy::getSimpleSphereIteratorPy,
        "Returns an iterator for enumerating spheres in this collection.\n"
        "@rtype: L{SimpleSphereCollectionIterator}\n"
        "@return: Sphere iterator.",
        boost::python::return_value_policy<
          boost::python::return_by_value,
          boost::python::with_custodian_and_ward_postcall<0, 1>
        >()
      )
      .def(
        "__len__",
        &SimpleSphereCollectionPy::getNumParticles,
        "Returns number of spheres in this collection.\n"
        "@rtype: int\n"
        "@return: number of spheres."
      )
      .def(
        "create",
        &SimpleSphereCollectionPy::createParticlePy,
        ( arg("sphere") ),
        boost::python::return_internal_reference<>(),
        "Creates a sphere in this collection, lifetime of sphere\n"
        "is tied to this collection.\n"
        "@type sphere: L{SimpleSphere}\n"
        "@kwarg sphere: A copy of this sphere is created.\n"
        "@rtype: L{SimpleSphere}\n"
        "@return: reference to newly created sphere."
      )
      .def(
        "translate",
        &SimpleSphereCollectionPy::translateByPy,
        ( arg("translation") ),
        "Translates all spheres in this collection by the specified amount.\n"
        "@type translation: L{esys.lsm.util.Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg translation: spheres are translated by this amount."
      )
      .def(
        "rotate",
        &SimpleSphereCollectionPy::rotatePy,
        ( arg("axis"), arg("axisPt") ),
        "Rotates all spheres in this collection about a specifed axis\n"
        "of rotation which passes through a specified point.\n"
        "@type axis: L{esys.lsm.util.FoundationPy.Vec3}\n"
        "@kwarg axis: axis of rotation, all spheres are rotated"
        " by C{axis.norm()}\n"
        " radians about the axis whose direction is C{axis/axis.norm()}.\n"
        "@type axisPt: L{esys.lsm.util.FoundationPy.Vec3}\n"
        "@kwarg axisPt: the axis of rotation is assumed to pass though"
        " this point."
      )
      .def_pickle(SimpleSphereCollectionPyPickleSuite())
      ;

      SimpleSphereCollectionPy::SimpleSphereIteratorPy::exportIterator(
        "SimpleSphereCollectionIterator"
      );
    }
  }
}
