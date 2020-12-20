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


#ifndef ESYS_LSMSTRESSTENSOR_H
#define ESYS_LSMSTRESSTENSOR_H

#include "Foundation/vec3.h"
#include "Foundation/Matrix3.h"
#include "Tools/StressCalculator/EigenvalueCalculator.h"
#include "Tools/StressCalculator/Contact.h"

#include <vector>
#include <map>

namespace esys
{
  namespace lsm
  {
    class Tensor
    {
    public:
      typedef EigenvalueCalculator::Complex       Complex;
      typedef EigenvalueCalculator::ComplexVector ComplexVector;

      Tensor()
        : m_pos(),
          m_tensor()
      {
      }

      Tensor(const Vec3 &pos, const Matrix3 &tensor)
        : m_pos(pos),
          m_tensor(tensor)
      {
      }

      virtual ~Tensor()
      {
      }

      const Vec3 &getPos() const
      {
        return m_pos;
      }

      const Matrix3 &getTensor() const
      {
        return m_tensor;
      }

      ComplexVector getEigenvalues() const
      {
        return EigenvalueCalculator().getEigenvalues(getTensor());
      }

    private:
      Vec3    m_pos;
      Matrix3 m_tensor;
    };

    class StressTensor : public Tensor
    {
    public:
      typedef EigenvalueCalculator::Complex       Complex;
      typedef EigenvalueCalculator::ComplexVector ComplexVector;

      StressTensor(const Vec3 &pos, const Matrix3 &tensor, double radius)
        : Tensor(pos, tensor),
          m_radius(radius)
      {
      }
      
      StressTensor(const ParticleData &particleData, const Matrix3 &tensor)
        : Tensor(particleData.getPos(), tensor),
          m_radius(particleData.getRad())
      {
      }
      
      virtual ~StressTensor()
      {
      }

      double getRad() const
      {
        return m_radius;
      }

    private:
      double m_radius;
    };
  }
}

#endif
