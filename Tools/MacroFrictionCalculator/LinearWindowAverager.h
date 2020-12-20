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


#ifndef ESYS_LSMLINEARWINDOWAVERAGER_H
#define ESYS_LSMLINEARWINDOWAVERAGER_H

#include <cstdlib>
#include "Foundation/vec3.h"

#include <vector>

namespace esys
{
  namespace lsm
  {
    class LinearWindowAverager
    {
    public:
      typedef std::vector<double> ValueVector;

      LinearWindowAverager(
        const ValueVector &valVector,
        int halfWindowSize,
        int beginIndex,
        int endIndex,
        int skipSize
      )
        :
          m_valVector(valVector),
          m_avValVector(),
          m_halfWindowSize(halfWindowSize),
          m_beginIndex(beginIndex),
          m_endIndex(min(endIndex, static_cast<int>(valVector.size()))),
          m_skipSize(skipSize)      
      {
      }

      const ValueVector &getAveragedVector()
      {
        if (m_avValVector.size() <= 0) {
          calculateAverageVals();
        }
        return m_avValVector;
      }

    protected:
      void calculateAverageVals()
      {
        m_avValVector.clear();
        const int halfWindowSize = max(0, m_halfWindowSize);
        for (int i = m_beginIndex; i < m_endIndex; i += m_skipSize)
        {
          double avSum     = 0.0;
          double weightSum = 0.0;
          const int endJ = min(static_cast<int>(m_valVector.size()), i + halfWindowSize+1);
          for (int j = max(0, i - halfWindowSize); j < endJ; j++)
          {
            const double weight = (1.0 - (static_cast<double>(abs(i-j))/(halfWindowSize+1)));
            weightSum += weight;
            avSum += weight*m_valVector[j];
          }
          m_avValVector.push_back(avSum/weightSum);
        }
      }

    private:
      ValueVector m_valVector;
      ValueVector m_avValVector;
      int         m_halfWindowSize;
      int         m_beginIndex;
      int         m_endIndex;
      int         m_skipSize;
    };
  }
}

#endif
