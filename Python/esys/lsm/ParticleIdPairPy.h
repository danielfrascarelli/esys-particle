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

#ifndef ESYS_LSMPARTICLEIDPAIRPY_H
#define ESYS_LSMPARTICLEIDPAIRPY_H

#include <boost/python.hpp>
#include <Python/esys/lsm/util/SetPy.h>

namespace esys
{
  namespace lsm
  {
    typedef std::pair<int,int> ParticleIdPair;

    class ParticleIdPairPy : public ParticleIdPair
    {
    public:
      typedef ParticleIdPair Inherited;
      
      ParticleIdPairPy(int id1, int id2);

      ParticleIdPairPy(const Inherited &pair);

      bool operator<(const ParticleIdPairPy &pair) const;

      int len() const;

      int getItem(int i);

      long hash() const;

      class PickleSuite : public boost::python::pickle_suite
      {
      public:
        static
        boost::python::tuple
        getinitargs(ParticleIdPairPy const& s);
      };

    };

    void exportParticleIdPair();
  }
}

#endif
