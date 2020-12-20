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


#ifndef ESYS_LSMEIGENVALUECALCULATOR_H
#define ESYS_LSMEIGENVALUECALCULATOR_H

#include "Foundation/Matrix3.h"

#include <iostream>
#include <complex>
#include <vector>
#include <algorithm>

namespace esys
{
  namespace lsm
  {
    class EigenvalueCalculator
    {
    public:
      typedef std::complex<double> Complex;
      typedef std::vector<Complex> ComplexVector;

      static const double ONE_THIRD;
      static const double SQRT_THREE;

      inline double cbrt(double a)
      {
          // cube root is antisymmetric
          return (a < 0) ? -std::pow(-a, ONE_THIRD) : std::pow(a, ONE_THIRD);
      }
      
      bool hasComplexEigenvalues(const Matrix3 &a) const
      {
          double b = -(a(0, 0) + a(1, 1) + a(2, 2));
          
          double a00a11 = a(0, 0) * a(1, 1);
          double a11a22 = a(1, 1) * a(2, 2);
          double a00a22 = a(0, 0) * a(2, 2);
          double a10a01 = a(1, 0) * a(0, 1);
          double a12a21 = a(1, 2) * a(2, 1);
          double a20a02 = a(2, 0) * a(0, 2);
          double c = a00a11 + a11a22 + a00a22 - a10a01 - a12a21 - a20a02;
      
          double a10a01a22 = a10a01 * a(2, 2);
          double a12a21a00 = a12a21 * a(0, 0);
          double a20a02a11 = a20a02 * a(1, 1);
          double a00a11a22 = a00a11 * a(2, 2);
          double a01a12a20 = a(0, 1) * a( 1, 2) * a( 2, 0);
          double a02a10a21 = a(0, 2) * a( 1, 0) * a( 2, 1);
          double d = a10a01a22 + a12a21a00 + a20a02a11
                   - a00a11a22 - a01a12a20 - a02a10a21;
      
          double bb = b * b;
          double Q = (3.0 * c - bb) / 9.0;
          double R = (9.0 * b * c - 27.0 * d - 2.0 * bb * b) / 54.0;
      
          return (Q * Q * Q > -(R * R));
      }
      
      class ComplexRealImagComparer
      {
      public:
        bool operator()(const Complex &c1, const Complex &c2) const
        {
          return 
            (
              (c1.real() < c2.real())
              ||
              (
                (c1.real() == c2.real())
                &&
                (c1.imag() < c2.imag())
              )
            );
        }
      };

      class ComplexAbsRealImagComparer
      {
      public:
        bool operator()(const Complex &c1, const Complex &c2) const
        {
          return 
            (
              (fabs(c1.real()) < fabs(c2.real()))
              ||
              (
                (fabs(c1.real()) == fabs(c2.real()))
                &&
                (fabs(c1.imag()) < fabs(c2.imag()))
              )
            );
        }
      };

      class ComplexNormComparer
      {
      public:
        bool operator()(const Complex &c1, const Complex &c2) const
        {
          return (std::norm(c1) < std::norm(c2));
        }
      };

      ComplexVector getEigenvalues(const Matrix3 &a)
      {
          ComplexVector e(3, Complex(0.0, 0.0));
          double b = -(a(0, 0) + a(1, 1) + a(2, 2));
          
          double a00a11 = a(0, 0) * a(1, 1);
          double a11a22 = a(1, 1) * a(2, 2);
          double a00a22 = a(0, 0) * a(2, 2);
          double a10a01 = a(1, 0) * a(0, 1);
          double a12a21 = a(1, 2) * a(2, 1);
          double a20a02 = a(2, 0) * a(0, 2);
          double c = a00a11 + a11a22 + a00a22 - a10a01 - a12a21 - a20a02;
      
          double a10a01a22 = a10a01 * a(2, 2);
          double a12a21a00 = a12a21 * a(0, 0);
          double a20a02a11 = a20a02 * a(1, 1);
          double a00a11a22 = a00a11 * a(2, 2);
          double a01a12a20 = a(0, 1) * a(1, 2) * a(2, 0);
          double a02a10a21 = a(0, 2) * a(1, 0) * a(2, 1);
          double d = a10a01a22 + a12a21a00 + a20a02a11
                   - a00a11a22 - a01a12a20 - a02a10a21;
      
          double bb = b * b;
          double Q = (3.0 * c - bb) / 9.0;
          double R = (9.0 * b * c - 27.0 * d - 2.0 * bb * b) / 54.0;
          double D = Q * Q * Q + R * R;
      
          double nOTb = -ONE_THIRD * b;
      
          if (D > 0) {
              // one real and two Complex conjugate
              double sqrtD = sqrt(D);
              double S = cbrt(R + sqrtD);
              double T = cbrt(R - sqrtD);
              double SpT = S + T;
              e[0] = nOTb + SpT;
      
              double x = nOTb - 0.5 * SpT;
              double y = 0.5 * SQRT_THREE * (S - T);
              e[1] = Complex(x, y);
              e[2] = Complex(x, -y);
          }
          else if (D < 0) {
              // all real and unequal
              double sqrtD = sqrt(-D);
              Complex S = std::pow(Complex(R, sqrtD), ONE_THIRD);  // T == S.conj()
              e[0] = nOTb + 2.0 * S.real();
      
              double x = nOTb - S.real();
              double y = SQRT_THREE * S.imag();
              e[1] = x - y;
              e[2] = x + y;
          }
          else {
              // all real and at least two equal
              double S = cbrt(R);  // T == S
              e[0] = nOTb + 2.0 * S;
              e[1] = nOTb - S;
              e[2] = e[1];
          }
          std::sort(e.begin(), e.end(), ComplexRealImagComparer());
          return e;
      }

      void printEigenvalues(const Matrix3 &a)
      {
          std::cout << "a: " << a(0, 0) << " " << a(0, 1) << " " << a(0, 2) << std::endl
                    << "   " << a(1, 0) << " " << a(1, 1) << " " << a(1, 2) << std::endl
                    << "   " << a(2, 0) << " " << a(2, 1) << " " << a(2, 2) << std::endl;
          
          if (hasComplexEigenvalues(a))
              std::cout << "Complex eigenvalues" << std::endl;
          else
              std::cout << "Real eigenvalues" << std::endl;

          ComplexVector e = getEigenvalues(a);
          std::cout << "e0: " << e[0] << std::endl;
          std::cout << "e1: " << e[1] << std::endl;
          std::cout << "e2: " << e[2] << std::endl;
          std::cout << std::endl;
      }
    };
  }
}

inline std::ostream &operator<<(std::ostream &oStream, const esys::lsm::EigenvalueCalculator::ComplexVector &vec)
{
  esys::lsm::EigenvalueCalculator::ComplexVector::const_iterator it = vec.begin();
  if (it != vec.end()) {
    oStream << (*it);
    it++;
  }
  for (; it != vec.end(); it++)
  {
    oStream << " " << (*it);
  }
  return oStream;
}

#endif
