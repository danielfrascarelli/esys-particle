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


#ifndef ESYS_LSMSTRESSTENSORCACLULATOR_H
#define ESYS_LSMSTRESSTENSORCACLULATOR_H

#include "Foundation/vec3.h"
#include "Foundation/Matrix3.h"
#include "Tools/StressCalculator/StressTensor.h"

#include <vector>
#include <map>

namespace esys
{
  namespace lsm
  {
    class IntercentreStressTensorCalculator
    {
    public:

      IntercentreStressTensorCalculator()
      {
      }

      template <typename TmplContactReference>
      void updateTensor(Matrix3 &tensor, const Vec3 &tensorPos, TmplContactReference &contact)
      {
          const Vec3 r = contact.getCentrePos2() - tensorPos;
          for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
              tensor(i,j) += contact.getForce()[i] * r[j];
            }
          }
      }
      
      template <typename TmplContactIterator>
      StressTensor calculate(TmplContactIterator &it)
      {
        Matrix3 tensor;
        ParticleData tensorData;
        if (it.hasNext())
        {
          tensorData = it.next().getParticle1();
          updateTensor(tensor, tensorData.getPos(), it.current());
        }
        
        while (it.hasNext())
        {
          updateTensor(tensor, tensorData.getPos(), it.next());
        }
        return StressTensor(tensorData, tensor);
      }
    };

    class ContactPtTensorCalculator
    {
    public:

      ContactPtTensorCalculator()
      {
      }

      template <typename TmplContactReference>
      void updateTensor(Matrix3 &tensor, const Vec3 &tensorPos, TmplContactReference &contact)
      {
        const Vec3 r = contact.getForcePos() - tensorPos;
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            tensor(i,j) += r[i] * contact.getForce()[j];
          }
        }
      }

      template <typename TmplContactIterator>
      StressTensor calculate(TmplContactIterator &it)
      {
        Matrix3 tensor;
        ParticleData tensorData;
        double volume = 1.0;
        if (it.hasNext())
        {
          tensorData = it.next().getParticle1();
          volume     = it.current().getVolume1();
          updateTensor(tensor, tensorData.getPos(), it.current());
        }
        while (it.hasNext())
        {
          updateTensor(tensor, tensorData.getPos(), it.next());
        }
        return StressTensor(tensorData, tensor/volume);
      }
    };    
  }
}

#endif
