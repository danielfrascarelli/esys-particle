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


#ifndef ESYS_LSMCLOSEPACKITERATOR_HPP
#define ESYS_LSMCLOSEPACKITERATOR_HPP

namespace esys
{
  namespace lsm
  {
    template <int NI, int NJ, int NK>
    TmplMatrix<NI,NJ,NK>::TmplMatrix()
    {
      for (int i = 0; i < NI; i++)
      {
        for (int j = 0; j < NJ; j++)
        {
          for (int k = 0; k < NJ; k++)
          {
            m_matrix[i][j][k] = 0.0;
          }
        }
      }
    }

    template <int NI, int NJ, int NK>
    TmplMatrix<NI,NJ,NK>::TmplMatrix(const TmplMatrix &m)
    {
      for (int i = 0; i < NI; i++)
      {
        for (int j = 0; j < NJ; j++)
        {
          for (int k = 0; k < NJ; k++)
          {
            m_matrix[i][j][k] = m(i,j,k);
          }
        }
      }
    }

    template <int NI, int NJ, int NK>
    TmplMatrix<NI,NJ,NK> &
    TmplMatrix<NI,NJ,NK>::TmplMatrix::operator=(const TmplMatrix &m)
    {
      for (int i = 0; i < NI; i++)
      {
        for (int j = 0; j < NJ; j++)
        {
          for (int k = 0; k < NJ; k++)
          {
            m_matrix[i][j][k] = m(i,j,k);
          }
        }
      }
      return *this;
    }

    template <int NI, int NJ, int NK>
    const double &TmplMatrix<NI,NJ,NK>::operator()(int i, int j, int k) const
    {
      return m_matrix[i][j][k];
    }

    template <int NI, int NJ, int NK>
    double &TmplMatrix<NI,NJ,NK>::operator()(int i, int j, int k)
    {
      return m_matrix[i][j][k];
    }

    template <int NI, int NJ, int NK>
    int TmplMatrix<NI,NJ,NK>::getNumI() const
    {
      return NI;
    }

    template <int NI, int NJ, int NK>
    int TmplMatrix<NI,NJ,NK>::getNumJ() const
    {
      return NJ;
    }

    template <int NI, int NJ, int NK>
    int TmplMatrix<NI,NJ,NK>::getNumK() const
    {
      return NK;
    }

    ClosePackIterator::ClosePackIterator()
      : m_radius(),
        m_minPt(),
        m_offsetMatrix(),
        m_dimRepeat(),
        m_dimCount(),
        m_dimIdx(),
        m_dim()
    {
    }

    ClosePackIterator::ClosePackIterator(
      int numI,
      int numJ,
      int numK,
      double sphereRadius,
      ClosePackOrientation orientation
    )
      : m_radius(sphereRadius),
        m_minPt(),
        m_offsetMatrix(),
        m_dimRepeat(),
        m_dimCount(numI, numJ, numK),
        m_dimIdx(),
        m_dim(s_orientationDimMap[orientation])
    {
      for (int i = 0; i < 3; i++)
      {
        if (m_dimCount[i] <= 0)
        {
          m_dimCount[i] = 1;
        }
      }
    }

    void ClosePackIterator::setDimRepeat(const Vec3L &dimRepeat)
    {
      m_dimRepeat = dimRepeat;
    }

    void ClosePackIterator::setOffsetMatrix(const OffsetMatrix &offsetMatrix)
    {
      m_offsetMatrix = offsetMatrix;
    }

    double ClosePackIterator::getRadius() const
    {
      return m_radius;
    }
    
    bool ClosePackIterator::hasNext() const
    {
       return (m_dimIdx[2] < m_dimCount[2]);
    }

    const Vec3 &ClosePackIterator::getMinPt() const
    {
      return m_minPt;
    }

    double ClosePackIterator::getOffset(int i) const
    {
      return
        m_offsetMatrix(
          i,
          (m_dimIdx[(i+1) % 3]) % m_dimRepeat[i],
          (m_dimIdx[(i+2) % 3]) % m_dimRepeat[i]
        );
    }

    void ClosePackIterator::incrementDimIndex()
    {
      (m_dimIdx[0])++;
      if (m_dimIdx[0] >= m_dimCount[0])
      {
        m_dimIdx[0] = 0;
        (m_dimIdx[1])++;
        if (m_dimIdx[1] >= m_dimCount[1])
        {
          m_dimIdx[1] = 0;
          (m_dimIdx[2])++;
        }
      }
    }
    
    Vec3 ClosePackIterator::next()
    {
      Vec3 centre;
      centre[m_dim[0]] = getOffset(0) + m_dimIdx[0]*(m_radius*2.0);
      centre[m_dim[1]] = getOffset(1) + m_dimIdx[1]*(m_radius*SQRT_3);
      centre[m_dim[2]] = getOffset(2) + m_dimIdx[2]*(m_radius*SQRT_8_OVER_3);
      incrementDimIndex();
      return (getMinPt() + centre);
    }
  }
}

#endif
