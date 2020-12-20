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
#include "Python/esys/lsm/geometry/TaggedIdConnectionPy.h"

namespace esys
{
  namespace lsm
  {
    TaggedIdConnectionPy::TaggedIdConnectionPy(Id id1, Id id2, Tag tag)
      : BasicInteraction(id1, id2, tag)
    {
    }

    TaggedIdConnectionPy::TaggedIdConnectionPy(const TaggedIdConnectionPy &idConn)
      : BasicInteraction(idConn)
    {
    }

    class TaggedIdConnectionPyPickleSuite : public boost::python::pickle_suite
    {
    public:
      static
      boost::python::tuple
      getinitargs(TaggedIdConnectionPy const& conn)
      {
          return boost::python::make_tuple(conn.getP1Id(), conn.getP2Id(), conn.getTag());
      }
    };

    using boost::python::arg;
    void exportTaggedIdConnection()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<TaggedIdConnectionPy>(
        "TaggedIdConnection",
        "Instances represent some kind of relationship (connection)\n"
        "between a pair of Id's.",
        boost::python::init<int, int, int>(
          ( arg("id1"), arg("id2"), arg("tag")=0 )
        )
      )
      .def(boost::python::init<const TaggedIdConnectionPy &>())
      .def("getId1", &TaggedIdConnectionPy::first)
      .def("getId2", &TaggedIdConnectionPy::second)
      .def("getTag", &TaggedIdConnectionPy::getTag)
      .def_pickle(TaggedIdConnectionPyPickleSuite())
      ;
    }
  }
}
