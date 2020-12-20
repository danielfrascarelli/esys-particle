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

#ifndef ESYS_LSMBONDINTERACTIONGROUPPY_H
#define ESYS_LSMBONDINTERACTIONGROUPPY_H

#include <Python/esys/lsm/InteractionGroupPy.h>
#include "Python/esys/lsm/ParticleIdPairSetPy.h"

namespace boost
{
  namespace python
  {
    class object;
  }
}

namespace esys
{
  namespace lsm
  {
    class LsmMpiPy;
    
    /**
     * Bond interaction group, delegates to an LsmMpiPy object.
     */
    class BondInteractionGroupPy : public InteractionGroupPy
    {
    public:

      BondInteractionGroupPy(
        LsmMpiPy          &lsm,
        const std::string &name
      );

      /**
       * Creates a bond between particles with specified ID's.
       * @param id1 Particle ID.
       * @param id2 Particle ID.
       */
      void createInteraction(int id1, int id2);

      /**
       * Creates bonds between specified pairs of particles.
       *
       * @param iterable Sequence of particle ID pairs.
       */
      void createInteractions(boost::python::object &iterable);

      /**
       * Return particle-id pairs indicating which particles
       * are bonded to one another.
       */
      ParticleIdPairSetPy getIdPairSet();
    };

    void exportBondInteractionGroup();

  } // namespace lsm
} // namespace esys

#endif // ESYS_LSMBONDINTERACTIONGROUPPY_H
