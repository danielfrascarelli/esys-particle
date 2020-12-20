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


#ifndef ESYS_LSMSTRESSTENSORCOLLECTION_H
#define ESYS_LSMSTRESSTENSORCOLLECTION_H

#include "Foundation/vec3.h"
#include "Foundation/Matrix3.h"
#include "Foundation/StlIterator.h"
#include "Tools/StressCalculator/StressTensor.h"

#include <vector>
#include <map>

namespace esys
{
  namespace lsm
  {
    template <typename TmplStressTensorCalculator>
    class StressTensorCollection
    {
    public:
      typedef TmplStressTensorCalculator          StressCalculator;
      typedef std::vector<StressTensor>           StressTensorVector;
      typedef ForwardIterator<StressTensorVector> StressTensorIterator;

      StressTensorCollection(StressCalculator &stressCalculator)
        : m_stressTensorVector(),
          m_pStressCalculator(&stressCalculator)
      {
      }

      template <typename TmplContactIterator>
      void addContactIterator(TmplContactIterator it)
      {
        m_stressTensorVector.push_back(m_pStressCalculator->calculate(it));
      }

      template <typename TmplContactIteratorIterator>
      void addContactIterators(TmplContactIteratorIterator it)
      {
        while (it.hasNext())
        {
          addContactIterator(it.next());
        }
      }
      
      StressTensorIterator getIterator()
      {
        return StressTensorIterator(m_stressTensorVector);
      }
      
      int size() const
      {
        return m_stressTensorVector.size();
      }

    private:
      StressTensorVector  m_stressTensorVector;
      StressCalculator    *m_pStressCalculator;
    };
  }
}

#endif
