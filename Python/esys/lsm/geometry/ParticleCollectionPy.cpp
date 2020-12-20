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
#include "Python/esys/lsm/geometry/ParticleCollectionPy.h"
#include "Python/esys/lsm/geometry/SimpleSpherePy.h"

namespace esys
{
  namespace lsm
  {
    ParticleCollectionPy::ParticleCollectionPy()
      : Inherited()
    {
    }

    ParticleCollectionPy::ParticleCollectionPy(
      const Inherited &particleCollection
    )
      : Inherited(particleCollection)
    {
    }

    ParticleCollectionPy::ParticleCollectionPy(
      const ParticleCollectionPy &particleCollection
    )
      : Inherited(particleCollection)
    {
    }

    ParticleCollectionPy::ParticleCollectionPy(
      ParticlePoolPtr particlePoolPtr
    )
      : Inherited(particlePoolPtr)
    {
    }

    BoundingBoxPy ParticleCollectionPy::getParticleBBoxPy() const
    {
      return
        BoundingBoxPy(
          ParticleCollection<SimpleSpherePy>::getParticleBBox()
        );
    }

    ParticleCollectionPy::ParticleIteratorPy
    ParticleCollectionPy::getParticleIteratorPy()
    {
      return
        ParticleIteratorPy(
          ParticleCollection<SimpleSpherePy>::getParticleIterator()
        );
    }

    SimpleSpherePy &ParticleCollectionPy::createParticlePy(const SimpleSpherePy &p)
    {
      return ParticleCollection<SimpleSpherePy>::createParticle(p);
    }

    void ParticleCollectionPy::rotatePy(const Vec3Py &rotation, const Vec3Py &pt)
    {
      rotate(rotation, pt);
    }

    void ParticleCollectionPy::translateByPy(const Vec3Py &translation)
    {
      translateBy(translation);
    }

    class ParticleCollectionPyPickleSuite : public boost::python::pickle_suite
    {
    public:
      static
      boost::python::tuple
      getstate(boost::python::object pcObj)
      {
        const ParticleCollectionPy &pc =
          boost::python::extract<const ParticleCollectionPy &>(pcObj);
        boost::python::list l;
        ParticleCollectionPy::ParticleConstIterator it = pc.getParticleIterator();
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

        ParticleCollectionPy &pc =
          boost::python::extract<ParticleCollectionPy &>(pcObj);
        boost::python::list l = boost::python::extract<boost::python::list>(state[1]);
        for (int i = 0; i < bpu::len(l); i++)
        {
          pc.createParticle(
            boost::python::extract<ParticleCollectionPy::Particle &>(l[i])
          );
        }
      }

      static bool getstate_manages_dict()
      {
        return true;
      }
    };

    using boost::python::arg;
    void exportParticleCollection()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<ParticleCollectionPy>(
        "ParticleCollection",
        "A collection of L{SimpleSphere} objects.",
        boost::python::init<>()
      )
      .def(
        "getParticleBBox",
        &ParticleCollectionPy::getParticleBBoxPy,
        "Returns the axis aligned bounding box of all particles in\n"
        "this collection.\n"
        "@rtype: L{esys.lsm.util.FoundationPy.BoundingBox}\n"
        "@return: tight bounding box for all particles in this collection."
      )
      .def(
        "getNumParticles",
        &ParticleCollectionPy::getNumParticles,
        "Returns the number of particles in this collection.\n"
        "@rtype: int\n"
        "@return: number of particles"
      )
      .def(
        "__iter__",
        &ParticleCollectionPy::getParticleIteratorPy,
        "Returns an iterator for enumerating particles in this collection.\n"
        "@rtype: ParticleCollectionParticleIterator\n"
        "@return: particle iterator.",
        boost::python::return_value_policy<
          boost::python::return_by_value,
          boost::python::with_custodian_and_ward_postcall<0, 1>
        >()
      )
      .def(
        "__len__",
        &ParticleCollectionPy::getNumParticles,
        "Returns number of particles in this collection.\n"
        "@rtype: int\n"
        "@return: number of particles."
      )
      .def(
        "createParticle",
        &ParticleCollectionPy::createParticlePy,
        ( arg("particle") ),
        boost::python::return_internal_reference<>(),
        "Creates a particle in this collection, lifetime of particle\n"
        "is tied to this collection.\n"
        "@type particle: L{SimpleSphere}\n"
        "@kwarg particle: A copy of this particle is created.\n"
        "@rtype: L{SimpleSphere}\n"
        "@return: reference to newly created particle."
      )
      .def(
        "translate",
        &ParticleCollectionPy::translateByPy,
        ( arg("translation") ),
        "Translates all particles in this collection by the specified amount.\n"
        "@type translation: L{esys.lsm.util.FoundationPy.Vec3}\n"
        "@kwarg translation: particles are translated by this amount."
      )
      .def(
        "rotate",
        &ParticleCollectionPy::rotatePy,
        ( arg("axis"), arg("pt") ),
        "Rotates all particles in this collection about a specifed axis\n"
        "of rotation which passes through a specified point.\n"
        "@type axis: L{esys.lsm.util.FoundationPy.Vec3}\n"
        "@kwarg axis: axis of rotation, all particles are rotated"
        " by C{axis.norm()}\n"
        " radians about the axis whose direction is C{axis/axis.norm()}.\n"
        "@type pt: L{esys.lsm.util.FoundationPy.Vec3}\n"
        "@kwarg pt: the axis of rotation is assumed to pass though"
        " this point."
      )
      .def_pickle(ParticleCollectionPyPickleSuite())
      ;

//      ParticleCollectionPy::ParticleIteratorPy::exportIterator(
//        "ParticleCollectionParticleIterator"
//      );
    }
  }
}
