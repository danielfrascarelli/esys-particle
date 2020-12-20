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
#include "Python/esys/lsm/ParticleIdPairPy.h"

namespace esys
{
  namespace lsm
  {
    ParticleIdPairPy::ParticleIdPairPy(int id1, int id2)
      : Inherited(std::min(id1, id2), std::max(id1,id2))
    {
    }

    ParticleIdPairPy::ParticleIdPairPy(const Inherited &pair)
      : Inherited(
          std::min(pair.first, pair.second),
          std::max(pair.first, pair.second)
        )
    {
    }

    bool ParticleIdPairPy::operator<(const ParticleIdPairPy &pair) const
    {
      bool result = first < pair.first;
      if (!result && (first == pair.first))
      {
        result = second < pair.second;
      }
      return result;
    }

    int ParticleIdPairPy::len() const
    {
      return 2;
    }

    int ParticleIdPairPy::getItem(int i)
    {
      const int origI = i;
      if (i < 0)
      {
          i += len();
      }
      if ((i >= len()) || i < 0)
      {
        std::stringstream msg;
        msg << "Index " << origI << " out of range.";

        PyErr_SetString(PyExc_IndexError, msg.str().c_str());
        boost::python::throw_error_already_set();
      }

      int item = INT_MIN;
      if (i == 0)
      {
        item = first;
      }
      else if (i == 1)
      {
        item = second;
      }
      return item;
    }

    long ParticleIdPairPy::hash() const
    {
      return first*second;
    }

    boost::python::tuple
    ParticleIdPairPy::PickleSuite::getinitargs(ParticleIdPairPy const& p)
    {
      return boost::python::make_tuple(p.first, p.second);
    }

    void exportParticleIdPair()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<ParticleIdPairPy>(
        "ParticleIdPair",
        "Pair of particle Id's.",
        boost::python::init<int,int>(
          (
            boost::python::arg("id1"),
            boost::python::arg("id2")
          ),
          "Construct pair with specified id's.\n"
          "@type id1: int\n"
          "@kwarg id1: Particle Id.\n"
          "@type id2: int\n"
          "@kwarg id2: Particle Id.\n"
        )
      )
      .def(
        "__getitem__",
        &ParticleIdPairPy::getItem,
        (boost::python::arg("idx")),
        "Returns id of indexed element.\n"
        "@type idx: int\n"
        "@kwarg idx: index value 0 or 1.\n"
      )
      .def(
        "__len__",
        &ParticleIdPairPy::len,
        "Returns 2.\n"
        "@rtype: int\n"
        "@return: length of pair (always 2).\n"
      )
      .def(boost::python::self < boost::python::self)
      .def(boost::python::self == boost::python::self)
      .def("__hash__", &ParticleIdPairPy::hash)
      .def_pickle(ParticleIdPairPy::PickleSuite())
      ;
    }
  }
}
