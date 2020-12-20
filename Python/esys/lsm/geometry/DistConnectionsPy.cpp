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
#include <patchlevel.h>
#include "Python/esys/lsm/geometry/DistConnectionsPy.h"
#include "Python/esys/lsm/geometry/SimpleBlockPy.h"
#include "Python/esys/lsm/geometry/SimpleSpherePy.h"
#include "Python/esys/lsm/geometry/TaggedIdConnectionPy.h"
#include "Python/esys/lsm/geometry/IteratorPy.h"
#include "Geometry/DistConnections.h"

namespace esys
{
  namespace lsm
  {
    template <typename TmplExtractType>
    class PythonIterIterator
    {
    public:
      PythonIterIterator(boost::python::object iter)
        : m_hasNext(true),
          m_next(),
          m_iter(iter)
      {
        update();
      }

      bool hasNext() const
      {
        return m_hasNext;
      }

      TmplExtractType next()
      {
        boost::python::object next = m_next;
        update();
        return boost::python::extract<TmplExtractType>(next);
      }

      void update()
      {
        try
        {
#if PY_VERSION_HEX >= 0x03000000
          m_next = m_iter.attr("__next__")();
#else
          m_next = m_iter.attr("next")();
#endif
        }
        catch (boost::python::error_already_set &e)
        {
          if (!PyErr_ExceptionMatches(PyExc_StopIteration))
          {
            throw;
          }
          PyErr_Clear();
          m_hasNext = false;
        }
      }

    private:
      bool                  m_hasNext;
      boost::python::object m_next;
      boost::python::object m_iter;
    };
    
    class DistConnectionsPy
      : public DistConnections<SimpleSpherePy, TaggedIdConnectionPy>
    {
    public:
      DistConnectionsPy(
        double maxDist=0.01,
        Tag defaultTag = 0,
        boost::python::object iteratable = boost::python::object()
      )
        : DistConnections<SimpleSpherePy, TaggedIdConnectionPy>(
            maxDist,
            defaultTag
          )
      {
        if (iteratable != boost::python::object())
        {
          addParticles(iteratable);
        }
      }

      void pythonObjectAddParticles(boost::python::object &iteratable)
      {
        create(
          PythonIterIterator<SimpleSpherePy &>(iteratable.attr("__iter__")())
        );
      }

      void pythonObjectAddParticlesWithTag(
        boost::python::object &iteratable,
        int connectionTag
      )
      {
        create(
          PythonIterIterator<SimpleSpherePy &>(iteratable.attr("__iter__")()),
          connectionTag
        );
      }

      void addParticles(boost::python::object &iteratable)
      {
        pythonObjectAddParticles(iteratable);
      }

      void addParticlesWithTag(
        boost::python::object &iteratable,
        int connectionTag
      )
      {
        pythonObjectAddParticlesWithTag(iteratable, connectionTag);
      }

      typedef
        IteratorPy<
          DistConnections<SimpleSpherePy,TaggedIdConnectionPy>::Iterator
        >
        IteratorPy2;
      IteratorPy2 getIterator() const
      {
        return
          IteratorPy2(
            DistConnections<
              SimpleSpherePy,
              TaggedIdConnectionPy
            >::getIterator()
          );
      }
    };

    using boost::python::arg;
    void exportDistConnections()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<DistConnectionsPy, boost::noncopyable>(
        "ConnectionFinder",
        boost::python::init<
          boost::python::optional<
            double,
            int,
            boost::python::object &
          >
        >(
          (
            arg("maxDist"),
            arg("bondTag"),
            arg("pList")
          ),
          "@type maxDist: float\n"
          "@kwarg maxDist: Threshhold distance governing which particles"
          " are connected.\n"
          "@type bondTag: int\n"
          "@kwarg bondTag: Default tag given to the connections.\n"
          "@type pList: iterable\n"
          "@kwarg pList: Iterable whose next method returns a SimpleSphere"
          " object."
        )
      )
        .def(
          "getNumParticles",
          &DistConnectionsPy::getNumParticles,
          "Returns the number of particles which have been added to this"
          " connection finder.\n"
          "@rtype: int\n"
          "@return: The number of particles added to this connection finder."
        )
        .def(
          "getNumConnections",
          &DistConnectionsPy::getNumConnections,
          "Returns the number of connections found by this connection finder.\n"
          "@rtype: int\n"
          "@return: number of L{TaggedIdConnection} objects in this collection."
        )
        .def("__iter__", &DistConnectionsPy::getIterator)
        .def("__len__", &DistConnectionsPy::getNumConnections)
        .def(
          "addParticles",
          &DistConnectionsPy::pythonObjectAddParticles,
          (
            arg("pList")
          ),
          "Adds SimpleSphere objects to this connection finder and determines"
          " any new connections.\n"
          "@type pList: iterable\n"
          "@kwarg pList: Iterator returning L{SimpleSphere} objects."
        )
        .def(
          "addParticles",
          &DistConnectionsPy::pythonObjectAddParticlesWithTag,
          (
            arg("pList"),
            arg("bondTag")
          ),
          "Adds SimpleSphere objects to this connection finder and determines"
          " any new connections.\n"
          //"@type pList: iterable\n"
          //"@kwarg pList: Iterator returning L{SimpleSphere} objects.\n"
          "@type bondTag: int\n"
          "@kwarg bondTag: New connections are assigned this tag."
        )
      ;

      boost::python::class_<DistConnectionsPy, boost::noncopyable>(
        "DistConnections",
        boost::python::init<
          boost::python::optional<
            double,
            int,
            boost::python::object &
          >
        >(
          (
            arg("maxDist"),
            arg("connTag"),
            arg("iterable")
          ),
          "@type maxDist: float\n"
          "@kwarg maxDist: Threshhold distance governing which particles"
          " are connected.\n"
          "@type connTag: int\n"
          "@kwarg connTag: Default tag given to the connections.\n"
          "@type iterable: iterable\n"
          "@kwarg iterable: Iterable whose next method returns a SimpleSphere"
          " object."
        )
      )
      .def(
        "getNumParticles",
        &DistConnectionsPy::getNumParticles,
        "Returns the number of particles which have been added to this"
        " connection finder.\n"
        "@rtype: int\n"
        "@return: The number of particles added to this connection finder."
      )
      .def(
        "getNumConnections",
        &DistConnectionsPy::getNumConnections,
        "Returns the number of connections found by this connection finder.\n"
        "@rtype: int\n"
        "@return: number of L{TaggedIdConnection} objects in this collection."
      )
      .def("__iter__", &DistConnectionsPy::getIterator)
      .def("__len__", &DistConnectionsPy::getNumConnections)
      .def(
        "addParticles",
        &DistConnectionsPy::pythonObjectAddParticles,
        (
          arg("iterable")
        ),
        "Adds SimpleSphere objects to this connection finder and determines"
        " any new connections.\n"
        "@type iterable: iterable\n"
        "@kwarg iterable: Iterator returning L{SimpleSphere} objects."
      )
      .def(
        "addParticles",
        &DistConnectionsPy::pythonObjectAddParticlesWithTag,
        (
          arg("iterable"),
          arg("connTag")
        ),
        "Adds SimpleSphere objects to this connection finder and determines"
        " any new connections.\n"
        //"@type iterable: iterable\n"
        //"@kwarg iterable: Iterator returning L{SimpleSphere} objects.\n"
        "@type connTag: int\n"
        "@kwarg connTag: New connections are assigned this tag."
      )
      ;

      DistConnectionsPy::IteratorPy2::exportIterator(
        "TaggedIdConnectionIterator"
      );
    }
  }
}
