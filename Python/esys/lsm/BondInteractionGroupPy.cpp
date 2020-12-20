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

#include <mpi.h>
#include <boost/version.hpp>
#include <boost/python.hpp>
#include "Python/esys/lsm/BondInteractionGroupPy.h"
#include "Python/esys/lsm/LsmMpiPy.h"
#include "Python/BoostPythonUtil/PythonIterIterator.h"

namespace esys
{
  namespace lsm
  {
    typedef LsmMpiPy::ParticleIdPair       ParticleIdPair;
    typedef LsmMpiPy::ParticleIdPairVector ParticleIdPairVector;
    BondInteractionGroupPy::BondInteractionGroupPy(
      LsmMpiPy          &lsm,
      const std::string &name
    )
      : InteractionGroupPy(lsm, name)
    {
    }

    void BondInteractionGroupPy::createInteractions(
      boost::python::object &iterable
    )
    {
      //
      // Convert the ID-pairs in iterable to
      // std::vector< std::pair<int,int> >
      //
      bpu::PythonIterIterator<boost::python::object> it(iterable);
      ParticleIdPairVector idPairVector;
      while (it.hasNext())
      {
        boost::python::object pair = it.next();
        idPairVector.push_back(
          ParticleIdPair(
            boost::python::extract<int>(pair[0]),
            boost::python::extract<int>(pair[1])
          )
        );
      }
      getLsm().createBonds(getName(), idPairVector);
    }

    void BondInteractionGroupPy::createInteraction(int id1, int id2)
    {
      boost::python::list pairList;
      pairList.append(boost::python::make_tuple(id1, id2));
      createInteractions(pairList);
    }

    ParticleIdPairSetPy BondInteractionGroupPy::getIdPairSet()
    {
      ParticleIdPairVector idPairVector =
        getLsm().getBondGroupIdPairs(getName());
      
      ParticleIdPairSetPy idPairSet;
      idPairSet.insert(idPairVector.begin(), idPairVector.end());
      return idPairSet;
    }

    void exportBondInteractionGroup()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<
        BondInteractionGroupPy,
        boost::python::bases<InteractionGroupPy>
      >(
        "BondInteractionGroup",
        "Base class for bonded interaction groups.",
        boost::python::no_init
      )
      .def(
        "createInteraction",
        &BondInteractionGroupPy::createInteraction,
        (
          boost::python::arg("id1"),
          boost::python::arg("id2")
        ),
        "Creates a bond between particles with specified ID's.\n"
        "@type id1: int\n"
        "@kwarg id1: Particle ID.\n"
        "@type id2: int\n"
        "@kwarg id2: Particle ID.\n"
      )
      .def(
        "createInteractions",
        &BondInteractionGroupPy::createInteractions,
        (boost::python::arg("idPairIterable")),
        "Creates bonds between particles with specified ID's.\n"
        "@type idPairIterable: iterable\n"
        "@kwarg idPairIterable: Supports C{iter(idPairIterable)} and"
        " each element in an iteration is a particle id pair.\n"
      )
      .def(
        "getIdPairSet",
        &BondInteractionGroupPy::getIdPairSet,
        "Returns pairs of particle-id's indicating pair-wise bonds.\n"
        "@rtype: L{ParticleIdPairSet}\n"
        "@return: Set of L{ParticleIdPair} objects.\n"
      )
      ;
    }
  }
}
