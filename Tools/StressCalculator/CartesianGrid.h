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


#ifndef ESYS_LSMCARTESIANGRID_H
#define ESYS_LSMCARTESIANGRID_H

#include "Foundation/vec3.h"
#include "Geometry/Vec3L.h"
#include "Foundation/BoundingBox.h"
#include "Foundation/StlIterator.h"

#include <vector>
#include <map>
#include <stdexcept>
#include <sstream>

#include <boost/pool/object_pool.hpp>
#include <boost/shared_ptr.hpp>

namespace esys
{
  namespace lsm
  {
    template <typename TmplValue>
    class CartesianGrid
    {
    public:
      typedef TmplValue         value_type;
      typedef value_type&       reference;
      typedef const value_type& const_reference;
      typedef value_type*       pointer;

      class Cell
      {
      public:
        class PosValuePair
        {
        public:
          PosValuePair(const Vec3 &pos, reference value)
            : m_pos(pos),
              m_pValue(&value)
          {
          }
          
          const Vec3 &getPos() const
          {
            return m_pos;
          }
          
          const_reference getValue() const
          {
            return *m_pValue;
          }

          reference getValue()
          {
            return *m_pValue;
          }

        private:
          Vec3    m_pos;
          pointer m_pValue;
        };

        typedef std::vector<PosValuePair>            PosValueVector;
        typedef ForwardIterator<PosValueVector>      Iterator;
        typedef ForwardConstIterator<PosValueVector> ConstIterator;
        
        Cell() : m_pos(), m_pointerVector(), m_pGrid(NULL)
        {
        }

        void insertRef(reference value)
        {
          m_pointerVector.push_back(PosValuePair(getPos(), value));
        }

        void insertRef(const Vec3 &pos, reference value)
        {
          m_pointerVector.push_back(PosValuePair(pos, value));
        }

        void setPos(const Vec3 &pos)
        {
          m_pos = pos;
        }

        const Vec3 &getPos() const
        {
          return m_pos;
        }

        CartesianGrid &getGrid()
        {
          return *m_pGrid;
        }

        const CartesianGrid &getGrid() const
        {
          return *m_pGrid;
        }

        void setGrid(CartesianGrid &grid)
        {
          m_pGrid = &grid;
        }

        BoundingBox getBox() const
        {
          return 
            BoundingBox(
              getPos() - getGrid().getGridSpacing()/2.0,
              getPos() + getGrid().getGridSpacing()/2.0
            );
        }

        /**
         * Return the number of elements in this cell.
         */
        int size() const
        {
          return m_pointerVector.size();
        }

        /**
         * Clears all elements from this cell.
         */
        void clear()
        {
          m_pointerVector.clear();
        }

        ConstIterator getIterator() const
        {
          return ConstIterator(m_pointerVector);
        }

        Iterator getIterator()
        {
          return Iterator(m_pointerVector);
        }

      private:
        Vec3           m_pos;
        PosValueVector m_pointerVector;
        CartesianGrid  *m_pGrid;
      };

      class VecIndexIterator
      {
      public:
        VecIndexIterator(const Vec3L &minIdx, const Vec3L &maxIdx)
          : m_minIndex(minIdx),
            m_maxIndex(maxIdx),
            m_index(minIdx)
        {
        }
        
        bool hasNext() const
        {
          return (m_index[2] <= m_maxIndex[2]);
        }
        
        const Vec3L &current()
        {
          return m_index;
        }

        void increment()
        {
          m_index[0]++;          
          if (m_index[0] > m_maxIndex[0]) {
            m_index[0] = m_minIndex[0];
            m_index[1]++;
            if (m_index[1] > m_maxIndex[1]) {
              m_index[1] = m_minIndex[1];
              m_index[2]++;
            }
          }        
        }

        Vec3L next()
        {
          const Vec3L index = current();
          increment();
          return index;
        }

      private:
        Vec3L         m_minIndex;
        Vec3L         m_maxIndex;
        Vec3L         m_index;
      };

      template <typename TmplGridPointer, typename TmplCellRef, typename TmplCell>
      class TCellIterator
      {
      public:
        typedef TmplCell    value_type;
        typedef TmplCellRef reference;
        TCellIterator(const Vec3L &minIdx, const Vec3L &maxIdx, TmplGridPointer pGrid)
          : m_idxIt(minIdx, maxIdx),
            m_pGrid(pGrid)
        {
        }

        bool hasNext() const
        {
          return m_idxIt.hasNext();
        }
        
        reference next()
        {
          return m_pGrid->getCell(m_pGrid->getScalarIndex(m_idxIt.next()));
        }
      private:
        VecIndexIterator m_idxIt;
        TmplGridPointer  m_pGrid;
      };

      typedef TCellIterator<CartesianGrid *, Cell &, Cell>             CellIterator;
      typedef TCellIterator<const CartesianGrid *, const Cell &, Cell> CellConstIterator;

      typedef std::vector<pointer> ValueVector;
      class ValueIterator : public ForwardIterator<ValueVector>
      {
      public:
        ValueIterator(ValueVector &v) : ForwardIterator<ValueVector>(v)
        {
        }

        typename CartesianGrid::reference next()
        {
          return *(ForwardIterator<ValueVector>::next());
        }
      };

      class ValueConstIterator : public ForwardConstIterator<ValueVector>
      {
      public:
        ValueConstIterator(const ValueVector &v) : ForwardConstIterator<ValueVector>(v)
        {
        }

        typename CartesianGrid::const_reference next()
        {
          return *(ForwardConstIterator<ValueVector>::next());
        }
      };

      CartesianGrid(const BoundingBox &bBox, double gridSpacing)
        : m_bBox(bBox),
          m_gridSpacing(gridSpacing),
          m_dimensions(),
          m_minIndex(),
          m_maxIndex(),
          m_cellVector(),
          m_valuePoolPtr(new Pool(1024)),
          m_valueVector()
      {
        initialise(bBox, gridSpacing);
      }

      const Vec3L &getDimensions() const
      {
        return m_dimensions;
      }

      const BoundingBox &getBBox() const
      {
        return m_bBox;
      }

      const Vec3 &getMinPt() const
      {
        return getBBox().getMinPt();
      }

      int getScalarIndex(int xIdx, int yIdx, int zIdx) const
      {
        return 
          xIdx*m_dimensions.Z()*m_dimensions.Y()
          +
          yIdx*m_dimensions.Z()
          +
          zIdx;
      }

      Vec3L getVecIndex(const Vec3 &pt) const
      {
        const Vec3 relPos = Vec3((pt - getMinPt())/m_gridSpacing);
        const Vec3L index = Vec3L(int(nearbyint(relPos.X())), int(nearbyint(relPos.Y())), int(nearbyint(relPos.Z())));
        return getMinVecIndex().max(getMaxVecIndex().min(index));
      }
      
      Vec3 getPos(const Vec3L &index) const
      {
        return 
          (
            getMinPt()
            +
            (Vec3(index.X(), index.Y(), index.Z())*m_gridSpacing)
          );
      }

      int getScalarIndex(const Vec3L &index) const
      {
        return getScalarIndex(index.X(), index.Y(), index.Z()); 
      }

      int getScalarIndex(const Vec3 &pt) const
      {
        return getScalarIndex(getVecIndex(pt));
      }

      const Vec3L &getMinVecIndex() const
      {
        return m_minIndex;
      }

      const Vec3L &getMaxVecIndex() const
      {
        return m_maxIndex;
      }

      void insert(const Vec3 &pos, const_reference data)
      {
        insertRef(pos, *(m_valuePoolPtr->construct(data)));
      }

      const Cell &getCell(int scalarIndex) const
      {
        return m_cellVector[scalarIndex];
      }
      
      Cell &getCell(int scalarIndex)
      {
        return m_cellVector[scalarIndex];
      }

      const Cell &getCell(const Vec3 &pos) const
      {
        return getCell(getScalarIndex(pos));
      }
      
      Cell &getCell(const Vec3 &pos)
      {
        return getCell(getScalarIndex(pos));
      }

      CellIterator getCellIterator(const Vec3 &pos, double radius)
      {
        return 
          CellIterator(
            getVecIndex(pos - radius),
            getVecIndex(pos + radius),
            this
          );
      }

      CellConstIterator getCellIterator(const Vec3 &pos, double radius) const
      {
        return 
          CellConstIterator(
            getVecIndex(pos - radius),
            getVecIndex(pos + radius),
            this
          );
      }

      CellIterator getCellIterator(const Vec3 &pos)
      {
        return getCellIterator(pos, 0.0);
      }

      CellConstIterator getCellIterator(const Vec3 &pos) const
      {
        return getCellIterator(pos, 0.0);
      }

      CellIterator getCellIterator()
      {
        return 
          CellIterator(
            getMinVecIndex(),
            getMaxVecIndex(),
            this
          );
      }

      CellConstIterator getCellIterator() const
      {
        return 
          CellConstIterator(
            getMinVecIndex(),
            getMaxVecIndex(),
            this
          );
      }

      ValueIterator getValueIterator()
      {
        return ValueIterator(m_valueVector);
      }

      ValueConstIterator getValueIterator() const
      {
        return ValueConstIterator(m_valueVector);
      }

      size_t size() const
      {
        return m_valueVector.size();
      }

      double getGridSpacing() const
      {
        return m_gridSpacing;
      }

      void clear()
      {
        CellIterator cellIt = getCellIterator();
        while (cellIt.hasNext())
        {
          cellIt.next().clear();
        }
        m_valuePoolPtr = PoolPtr(new Pool(1024));
        m_valueVector.clear();
      }

    protected:
      void insertRef(const Vec3 &pos, reference data)
      {
        getCell(pos).insertRef(pos, data);
        m_valueVector.push_back(&data);
      }

      void initialise(const BoundingBox &bBox, double gridSpacing)
      {
        m_bBox        = bBox;
        m_gridSpacing = gridSpacing;

        const Vec3 dims  = m_bBox.getSizes()/gridSpacing;
        m_dimensions = Vec3L(int(nearbyint(dims[0])), int(nearbyint(dims[1])), int(nearbyint(dims[2])));
        m_dimensions = m_dimensions.max(Vec3L(1, 1, 1));

        m_cellVector.clear();
        m_cellVector.insert(
          m_cellVector.end(),
          m_dimensions.X()*m_dimensions.Y()*m_dimensions.Z(),
          Cell()
        );

        m_valueVector.clear();
        m_valueVector.reserve(m_cellVector.size());

        m_minIndex = Vec3L(0, 0, 0);
        m_maxIndex = (m_dimensions - 1);

        for (int i = getMinVecIndex()[0]; i <= getMaxVecIndex()[0]; i++) {
          for (int j = getMinVecIndex()[1]; j <= getMaxVecIndex()[1]; j++) {
            for (int k = getMinVecIndex()[2]; k <= getMaxVecIndex()[2]; k++) {
              const Vec3L index(i, j, k);
              const Vec3 pos = getPos(index);
              Cell &cell = getCell(getScalarIndex(index));
              if ((&(getCell(pos))) == (&(cell)))
              {
                cell.setPos(pos);
                cell.setGrid(*this);
              }
              else
              {
                std::stringstream sStream;
                sStream 
                  << "Couldn't set grid location for cell ("
                  << i << "," << j << "," << k << "), pos = " << pos;
                throw std::runtime_error(sStream.str());
              }
            }
          }
        }
      }

    private:
      typedef std::vector<Cell>              CellVector;
      typedef boost::object_pool<value_type> Pool;
      typedef boost::shared_ptr<Pool>        PoolPtr;

      BoundingBox m_bBox;
      double      m_gridSpacing;
      Vec3L       m_dimensions;
      Vec3L       m_minIndex;
      Vec3L       m_maxIndex;
      CellVector  m_cellVector;      
      PoolPtr     m_valuePoolPtr;
      ValueVector m_valueVector;
    };
  }
}

#endif
