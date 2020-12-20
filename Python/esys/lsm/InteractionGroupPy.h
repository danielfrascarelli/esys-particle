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

#ifndef ESYS_LSMINTERACTIONGROUPPY_H
#define ESYS_LSMINTERACTIONGROUPPY_H

#include <string>

namespace esys
{
  namespace lsm
  {
    class LsmMpiPy;
    
    /**
     * Base class for python-exposed interaction group.
     * Delegates calls to LsmMpiPy object.
     */
    class InteractionGroupPy
    {
    public:
      InteractionGroupPy(
        LsmMpiPy          &lsmMpi,
        const std::string &name
      );

      const std::string &getName() const
      {
        return m_name;
      }

    protected:
      LsmMpiPy &getLsm()
      {
        return *m_pLsm;
      }

      const LsmMpiPy &getLsm() const
      {
        return *m_pLsm;
      }

    private:
      LsmMpiPy    *m_pLsm;
      std::string m_name;
    };

    void exportInteractionGroup();

  } // namespace lsm
} // namespace esys

#endif // ESYS_LSMINTERACTIONGROUPPY_H
