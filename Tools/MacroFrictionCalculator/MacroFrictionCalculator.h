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


#ifndef ESYS_LSMMACROFRICTIONCALCULATOR_H
#define ESYS_LSMMACROFRICTIONCALCULATOR_H

#include "Foundation/vec3.h"

#include <vector>

namespace esys
{
  namespace lsm
  {
    class MacroFrictionCalculator
    {
    public:
      typedef std::pair<Vec3,Vec3> WallForcePair;
      typedef std::vector<double>  FrictionVector;

      MacroFrictionCalculator(int normalDimIndex, int shearDimIndex)
        : m_normalDimIndex(normalDimIndex),
          m_shearDimIndex(shearDimIndex),
          m_frictionVector()
      {
      }

      double getFriction(const WallForcePair &forcePair) const
      {
        const double normalStress = (forcePair.first[m_normalDimIndex] - forcePair.second[m_normalDimIndex]);
        if (normalStress != 0.0) {
          const double shearStress = forcePair.first[m_shearDimIndex] - forcePair.second[m_shearDimIndex];
          return shearStress/normalStress;
        }
        return 0.0;
      }

      void add(const WallForcePair &forcePair)
      {
        m_frictionVector.push_back(getFriction(forcePair));
      }

      template<typename TmplIterator>
      void add(TmplIterator it)
      {
        while (it.hasNext())
        {
          add(it.next());
        }
      }

      const FrictionVector &getFrictionVector() const
      {
        return m_frictionVector;
      }

    private:
      int            m_normalDimIndex;
      int            m_shearDimIndex;
      FrictionVector m_frictionVector;
    };
  }
}

#endif
