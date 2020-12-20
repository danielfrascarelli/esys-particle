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

#ifndef ESYS_LSMGAUSSIANGRIDDER_H
#define ESYS_LSMGAUSSIANGRIDDER_H

#include "Foundation/vec3.h"

#include <math.h>

namespace esys
{
  namespace lsm
  {
    class GaussianGridder
    {
    public:
      GaussianGridder(double stdDeviation)
        : m_stdDeviation(stdDeviation)
      {
      }

      static const double SQRT_2PI;
      
      double getWeight(const Vec3 &regPos, const Vec3 &irrPos) const
      {
        return exp( - ((irrPos-regPos).norm2())/(2.0*m_stdDeviation*m_stdDeviation))/(SQRT_2PI*m_stdDeviation);
      }

      template <typename TmplCartesianGrid>
      void transform(const TmplCartesianGrid &irregular, TmplCartesianGrid &regular) const
      {
        typename TmplCartesianGrid::CellIterator regIt = regular.getCellIterator();
        const double radius = 3.99*m_stdDeviation;
        while (regIt.hasNext())
        {
          const typename TmplCartesianGrid::Cell &regCell = regIt.next();
          typename TmplCartesianGrid::CellConstIterator irrIt = irregular.getCellIterator(regCell.getPos(), radius);
          typename TmplCartesianGrid::value_type value;
          value *= 0.0;
          double weightSum = 0.0;
          while (irrIt.hasNext())
          {
            typename TmplCartesianGrid::Cell::ConstIterator pairIt = irrIt.next().getIterator();
            while (pairIt.hasNext())
            {
              const typename TmplCartesianGrid::Cell::PosValuePair &pair = pairIt.next();
              const double weight = getWeight(regCell.getPos(), pair.getPos());
              weightSum += weight;
              value += (pair.getValue()*weight);
            }
          }
          if (weightSum != 0) {
            value /= weightSum;
          }
          regular.insert(regCell.getPos(), value);
        }
      }

    private:
      double m_stdDeviation;
    };
  }
}

#endif
